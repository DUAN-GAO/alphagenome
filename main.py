import argparse
import os
import sys
import json
import logging
from datetime import datetime
import myvariant
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 从环境变量读取 API Key
API_KEY = os.environ.get("ALPHAGENOME_API_KEY")
if not API_KEY:
    logger.error("未设置环境变量 ALPHAGENOME_API_KEY，请设置后再运行")
    sys.exit(1)


def rsid_to_variant_info(rsid, genome_build="hg19"):
    """输入 rsID，返回 variant 所需信息 (chrom, pos, ref, alt)"""
    mv = myvariant.MyVariantInfo()
    try:
        out = mv.getvariant(rsid, fields="dbsnp")
    except Exception as e:
        raise RuntimeError(f"查询 {rsid} 失败: {e}")

    if not out or "dbsnp" not in out:
        raise ValueError(f"{rsid} 未找到相关 dbSNP 信息")

    dbsnp = out["dbsnp"]
    vartype = dbsnp.get("vartype", "").lower()
    if vartype != "snv":
        raise ValueError(f"{rsid} 不是单核苷酸变异 (SNV)，实际类型: {vartype}")

    ref_base = dbsnp.get("ref")
    alt_base = dbsnp.get("alt")
    chrom = dbsnp.get("chrom")
    pos = dbsnp.get(genome_build, {}).get("start")

    if None in [ref_base, alt_base, chrom, pos]:
        missing = []
        if ref_base is None: missing.append("ref")
        if alt_base is None: missing.append("alt")
        if chrom is None: missing.append("chrom")
        if pos is None: missing.append(f"pos({genome_build})")
        raise ValueError(f"{rsid} 坐标或碱基信息不完整，缺失: {', '.join(missing)}")

    return f"chr{chrom}", pos, ref_base, alt_base


def parse_variant_line(line):
    """
    解析一行变异文本，格式如:
    chrX:40074676__G > A
    返回 (chrom, pos, ref, alt, name)
    """
    line = line.strip()
    if not line or line.startswith('#'):
        return None

    # 按双下划线分割
    if '__' not in line:
        raise ValueError(f"行格式错误，缺少 '__' 分隔符: {line}")
    left, right = line.split('__', 1)
    left = left.strip()
    right = right.strip()

    # 解析左边：chrom:pos
    if ':' not in left:
        raise ValueError(f"左边格式错误，缺少 ':' : {left}")
    chrom_part, pos_str = left.split(':', 1)
    chrom = chrom_part.strip()
    # 如果染色体没有"chr"前缀，自动添加（保持与Variant构造函数一致）
    if not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
    try:
        pos = int(pos_str.strip())
    except ValueError:
        raise ValueError(f"位置解析错误: {pos_str}")

    # 解析右边：ref > alt（可能有空格）
    if ' > ' not in right:
        raise ValueError(f"右边格式错误，缺少 ' > ' : {right}")
    ref, alt = right.split(' > ', 1)
    ref = ref.strip()
    alt = alt.strip()

    # 生成名称（用于输出文件名）
    name = f"{chrom}:{pos}_{ref}>{alt}"
    return chrom, pos, ref, alt, name


def score_variant_from_coords(chrom, pos, ref, alt, name, dna_model):
    """
    直接使用坐标进行打分，返回结果字典
    """
    variant = genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases=alt,
        name=name
    )

    # 定义预测区间（长度2048）
    sequence_length = 2048
    interval = variant.reference_interval.resize(sequence_length)

    # 定义 scorer
    scorer = variant_scorers.CenterMaskScorer(
        width=None,
        aggregation_type=variant_scorers.AggregationType.DIFF_SUM_LOG2,
        requested_output=dna_client.OutputType.RNA_SEQ,
    )

    # 打分
    try:
        score_result = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=[scorer],
            organism=dna_client.Organism.HOMO_SAPIENS,
        )
    except Exception as e:
        raise RuntimeError(f"模型打分失败: {e}")

    # 提取结果（假设 score_result[0].var 是得分）
    result = {
        "name": name,
        "chromosome": chrom,
        "position": pos,
        "ref": ref,
        "alt": alt,
        "score_var": score_result[0].var if hasattr(score_result[0], 'var') else None,
        "score_raw": str(score_result[0]) if not hasattr(score_result[0], 'var') else None,
        "timestamp": datetime.now().isoformat()
    }
    return result


def score_variant_by_rsid(rsid, dna_model):
    """通过 rsID 进行打分"""
    try:
        chrom, pos, ref, alt = rsid_to_variant_info(rsid)
        name = rsid
    except Exception as e:
        logger.error(f"转换 rsID 失败 ({rsid}): {e}")
        return None

    try:
        result = score_variant_from_coords(chrom, pos, ref, alt, name, dna_model)
        result["rsid"] = rsid  # 额外保存rsID
        return result
    except Exception as e:
        logger.error(f"打分失败 ({rsid}): {e}")
        return None


def process_variants(variants, dna_model, outdir):
    """
    通用处理函数：对每个变异进行打分并保存
    variants: 列表，每个元素是 (chrom, pos, ref, alt, name) 或 (name, chrom, pos, ref, alt) 的元组
    """
    os.makedirs(outdir, exist_ok=True)
    all_results = []

    for item in variants:
        if len(item) == 5:
            chrom, pos, ref, alt, name = item
        else:
            continue

        logger.info(f"处理变异 {name} ...")
        try:
            result = score_variant_from_coords(chrom, pos, ref, alt, name, dna_model)
            all_results.append(result)

            # 每个变异单独保存 JSON 文件
            out_file = os.path.join(outdir, f"{name}.json")
            with open(out_file, "w") as f:
                json.dump(result, f, indent=2)
            logger.info(f"已保存结果至 {out_file}")
        except Exception as e:
            logger.error(f"处理变异 {name} 失败: {e}")

    # 汇总保存所有结果
    if all_results:
        summary_file = os.path.join(outdir, "all_results.json")
        with open(summary_file, "w") as f:
            json.dump(all_results, f, indent=2)
        logger.info(f"汇总结果保存至 {summary_file}")
    else:
        logger.warning("没有成功处理任何变异")


def main_from_file(file_path, outdir):
    """从文件读取变异列表并处理"""
    logger.info(f"读取变异文件: {file_path}")
    variants = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            try:
                parsed = parse_variant_line(line)
                if parsed is None:
                    continue
                variants.append(parsed)
            except Exception as e:
                logger.warning(f"第 {line_num} 行解析失败: {line.strip()} - {e}")

    if not variants:
        logger.error("未找到有效的变异行，退出")
        return

    logger.info(f"成功解析 {len(variants)} 个变异")
    dna_model = dna_client.create(API_KEY)
    process_variants(variants, dna_model, outdir)


def main_from_rsids(rsids, outdir):
    """处理 rsID 列表"""
    dna_model = dna_client.create(API_KEY)
    variants = []
    for rsid in rsids:
        try:
            chrom, pos, ref, alt = rsid_to_variant_info(rsid)
            variants.append((chrom, pos, ref, alt, rsid))
        except Exception as e:
            logger.error(f"转换 rsID 失败 ({rsid}): {e}")
    if variants:
        process_variants(variants, dna_model, outdir)
    else:
        logger.error("没有有效的 rsID")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="使用 AlphaGenome 对变异进行打分")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--rsid", action="append",
                       help="dbSNP rsID，可多次使用或逗号分隔，如 --rsid rs123 --rsid rs456")
    group.add_argument("--input", type=str,
                       help="包含变异的文本文件，每行格式: chrX:40074676__G > A")
    parser.add_argument("--outdir", default=".",
                        help="输出目录（默认当前目录）")
    args = parser.parse_args()

    if args.input:
        main_from_file(args.input, args.outdir)
    elif args.rsid:
        # 处理 --rsid 参数：可能包含逗号分隔的多个
        rsids = []
        for item in args.rsid:
            for part in item.split(','):
                part = part.strip()
                if part:
                    rsids.append(part)
        main_from_rsids(rsids, args.outdir)
    else:
        parser.error("请指定 --rsid 或 --input")