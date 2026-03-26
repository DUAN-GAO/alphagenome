import argparse
import pandas as pd
import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import variant_scorers

# 固定 API Key（保持不变）
API_KEY = "AIzaSyD5Kht8QzCPkHeJ456_Tf_eBWirtKhmaRU"


def parse_variant_line(line):
    """
    解析一行变异文本，格式如:
    chrX:40074676__G > A
    返回 (chrom, pos, ref, alt, name)
    """
    line = line.strip()
    if not line or line.startswith('#'):
        return None

    if '__' not in line:
        raise ValueError(f"行格式错误，缺少 '__' 分隔符: {line}")
    left, right = line.split('__', 1)
    left = left.strip()
    right = right.strip()

    if ':' not in left:
        raise ValueError(f"左边格式错误，缺少 ':' : {left}")
    chrom_part, pos_str = left.split(':', 1)
    chrom = chrom_part.strip()
    if not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
    try:
        pos = int(pos_str.strip())
    except ValueError:
        raise ValueError(f"位置解析错误: {pos_str}")

    if ' > ' not in right:
        raise ValueError(f"右边格式错误，缺少 ' > ' : {right}")
    ref, alt = right.split(' > ', 1)
    ref = ref.strip()
    alt = alt.strip()

    name = f"{chrom}:{pos}_{ref}>{alt}"
    return chrom, pos, ref, alt, name


def extract_scalar_score(score_result):
    """
    从 score_variant 返回的结果中提取一个标量分数。
    处理 DataFrame 类型，默认取所有轨道的 nonzero_mean 的平均值。
    """
    if not score_result:
        raise ValueError("score_result 为空")
    res = score_result[0]

    # 如果 .var 已经是数值
    if hasattr(res, 'var') and isinstance(res.var, (int, float, np.number)):
        return res.var

    # 如果 .var 是 pandas DataFrame
    if hasattr(res, 'var') and isinstance(res.var, pd.DataFrame):
        df = res.var
        # 优先使用 nonzero_mean 列的平均值（该列在错误输出中可见）
        if 'nonzero_mean' in df.columns:
            return df['nonzero_mean'].mean()
        # 否则，对所有数值列取全局均值
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            return df[numeric_cols].mean().mean()
        # 最后兜底
        raise TypeError("无法从 DataFrame 中提取数值")
    
    # 如果 .var 是 list/array 且长度为1，取第一个元素
    if hasattr(res, 'var') and hasattr(res.var, '__len__') and len(res.var) == 1:
        return res.var[0]
    
    # 如果 .score 存在且为数值
    if hasattr(res, 'score') and isinstance(res.score, (int, float)):
        return res.score

    raise RuntimeError(f"无法从结果中提取标量分数，结果类型: {type(res)}")


def main(input_file, outdir="."):
    print('API reached...')
    dna_model = dna_client.create(API_KEY)

    variants = []
    with open(input_file, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            try:
                parsed = parse_variant_line(line)
                if parsed is None:
                    continue
                variants.append(parsed)
            except Exception as e:
                print(f"警告：第 {line_num} 行解析失败: {line.strip()} - {e}")

    if not variants:
        print("错误：未找到有效的变异行，退出")
        return

    print(f"成功解析 {len(variants)} 个变异，开始打分...")

    all_results = []

    for chrom, pos, ref, alt, name in variants:
        print(f"处理变异 {name} ...")
        try:
            variant = genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt,
                name=name
            )

            sequence_length = 16384
            interval = variant.reference_interval.resize(sequence_length)

            scorer = variant_scorers.CenterMaskScorer(
                width=None,
                aggregation_type=variant_scorers.AggregationType.DIFF_SUM_LOG2,
                requested_output=dna_client.OutputType.RNA_SEQ,
            )

            score_result = dna_model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=[scorer],
                organism=dna_client.Organism.HOMO_SAPIENS,
            )

            score_value = extract_scalar_score(score_result)
            all_results.append({
                "name": name,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "score": score_value
            })
            print(f"  得分: {score_value}")
        except Exception as e:
            print(f"  错误: {e}")
            all_results.append({
                "name": name,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "score": f"ERROR: {e}"
            })

    output_file = f"{outdir}/results.txt" if outdir != "." else "results.txt"
    import os
    os.makedirs(outdir, exist_ok=True)
    with open(output_file, "w") as f:
        f.write("name\tchromosome\tposition\tref\talt\tscore\n")
        for r in all_results:
            f.write(f"{r['name']}\t{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t{r['score']}\n")

    print(f"[RESULT] 所有结果已保存到 {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="从文本文件读取变异并打分（AlphaGenome）")
    parser.add_argument("--input", type=str, required=True,
                        help="包含变异的文本文件，每行格式: chrX:40074676__G > A")
    parser.add_argument("--outdir", default=".",
                        help="输出目录（默认当前目录），结果文件名为 results.txt")
    args = parser.parse_args()
    main(args.input, args.outdir)