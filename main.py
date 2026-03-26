import argparse
import pandas as pd
import numpy as np
import os
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers

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


def extract_scalar_score(score_result, variant_name, agg_col='nonzero_mean', save_raw=False, outdir="."):
    """
    从 score_variant 返回的结果中提取标量分数。
    agg_col: 用于聚合的列名，默认为 'nonzero_mean'。
    """
    if not score_result:
        raise ValueError("score_result 为空")
    res = score_result[0]

    # 如果是数值，直接返回
    if hasattr(res, 'var') and isinstance(res.var, (int, float, np.number)):
        return res.var

    # 如果是 DataFrame
    if hasattr(res, 'var') and isinstance(res.var, pd.DataFrame):
        df = res.var
        if save_raw:
            raw_file = os.path.join(outdir, f"{variant_name}_raw.csv")
            df.to_csv(raw_file, index=False)
            print(f"    原始数据已保存至 {raw_file}")

        # 尝试使用指定列
        if agg_col in df.columns:
            score = df[agg_col].mean()
            # 打印该列的基本统计信息以帮助调试
            print(f"    {agg_col} 范围: {df[agg_col].min():.4f} - {df[agg_col].max():.4f}, 均值: {score:.4f}")
            return score
        else:
            # 如果指定列不存在，列出可用列并尝试使用第一个数值列
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) == 0:
                raise TypeError("DataFrame 中没有数值列")
            print(f"    警告: 列 '{agg_col}' 不存在，可用数值列: {list(numeric_cols)}")
            col = numeric_cols[0]
            score = df[col].mean()
            print(f"    使用列 '{col}' 的平均值: {score:.4f}")
            return score

    # 如果是数组
    if hasattr(res, 'var') and hasattr(res.var, '__len__') and len(res.var) == 1:
        return res.var[0]

    if hasattr(res, 'score') and isinstance(res.score, (int, float)):
        return res.score

    raise RuntimeError(f"无法从结果中提取标量分数，结果类型: {type(res)}")


def main(input_file, outdir=".", agg_col='nonzero_mean'):
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
    print(f"聚合列: {agg_col}")

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

            score_value = extract_scalar_score(score_result, name, agg_col=agg_col, save_raw=True, outdir=outdir)
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
    parser.add_argument("--agg_col", default="nonzero_mean",
                        help="用于聚合的列名（在 DataFrame 中），默认为 'nonzero_mean'")
    args = parser.parse_args()
    main(args.input, args.outdir, args.agg_col)