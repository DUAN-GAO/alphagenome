import argparse
import pandas as pd
import numpy as np
import os
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers

# 固定 API Key
API_KEY = "AIzaSyD5Kht8QzCPkHeJ456_Tf_eBWirtKhmaRU"


def parse_variant_line(line):
    """解析一行变异文本，格式如: chrX:40074676__G > A"""
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


def extract_scalar_score(score_result, variant_name, save_raw=False, outdir="."):
    """
    从 score_variant 返回的结果中提取标量分数。
    如果 save_raw 为 True，则将原始的 var DataFrame 保存到文件。
    """
    if not score_result:
        raise ValueError("score_result 为空")
    res = score_result[0]

    # 1) 如果是数值，直接返回
    if hasattr(res, 'var') and isinstance(res.var, (int, float, np.number)):
        return res.var

    # 2) 如果是 DataFrame
    if hasattr(res, 'var') and isinstance(res.var, pd.DataFrame):
        df = res.var
        # 可选：保存原始 DataFrame 以便检查
        if save_raw:
            raw_file = os.path.join(outdir, f"{variant_name}_raw.csv")
            df.to_csv(raw_file, index=False)
            print(f"    原始数据已保存至 {raw_file}")
        
        # 使用 nonzero_mean 列的平均值作为最终得分
        if 'nonzero_mean' in df.columns:
            score = df['nonzero_mean'].mean()
            return score
        else:
            # 如果不存在 nonzero_mean，则取所有数值列的平均值
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) > 0:
                score = df[numeric_cols].mean().mean()
                return score
            else:
                raise TypeError("DataFrame 中没有数值列")

    # 3) 其他情况（如数组）
    if hasattr(res, 'var') and hasattr(res.var, '__len__') and len(res.var) == 1:
        return res.var[0]

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

    os.makedirs(outdir, exist_ok=True)
    all_results = []

    for chrom, pos, ref, alt, name in variants:
        print(f"处理变异 {name} ...")
        print(f"  坐标: {chrom}:{pos}  {ref}>{alt}")

        try:
            variant = genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt,
                name=name
            )

            # 确认区间
            sequence_length = 16384
            interval = variant.reference_interval.resize(sequence_length)
            print(f"  区间: {interval.start} - {interval.end}")

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

            # 提取标量分数，并可选保存原始 DataFrame
            score_value = extract_scalar_score(score_result, name, save_raw=True, outdir=outdir)
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

    # 写入汇总结果
    output_file = os.path.join(outdir, "results.txt")
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