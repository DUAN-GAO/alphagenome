#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import myvariant
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import variant_scorers


# 固定 API Key
API_KEY = "AIzaSyD5Kht8QzCPkHeJ456_Tf_eBWirtKhmaRU"


def load_rsids_from_vcf(vcf_path):
    """从 VCF 文件提取 rsID"""
    opener = gzip.open if vcf_path.endswith(".gz") else open

    rsids = []
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            rsid = parts[2]
            if rsid.startswith("rs"):
                rsids.append(rsid)

    return rsids


def rsid_to_variant_info(rsid, genome_build="hg19"):
    """
    输入 rsID，返回:
    chromosome='chr22', position=36201698, reference_bases='A', alternate_bases='C'
    仅适用 SNV
    """
    mv = myvariant.MyVariantInfo()
    out = mv.getvariant(rsid, fields="dbsnp")
    if not out or "dbsnp" not in out:
        raise ValueError(f"{rsid} 未找到相关 dbSNP 信息")

    dbsnp = out["dbsnp"]
    vartype = dbsnp.get("vartype", "").lower()
    if vartype != "snv":
        raise ValueError(f"{rsid} 不是单核苷酸变异 (SNV)")

    ref_base = dbsnp.get("ref")
    alt_base = dbsnp.get("alt")
    chrom = dbsnp.get("chrom")
    pos = dbsnp.get(genome_build, {}).get("start")

    if None in [ref_base, alt_base, chrom, pos]:
        raise ValueError(f"{rsid} 坐标或碱基信息不完整")

    return f"chr{chrom}", pos, ref_base, alt_base


def score_single_rsid(rsid, dna_model):
    """使用 AlphaGenome 对单个 rsID 打分"""

    print(f"\n[INFO] 处理 {rsid}")

    # === rsID → chr/pos/ref/alt ===
    chrom, pos, ref, alt = rsid_to_variant_info(rsid)
    print(f"[INFO] {rsid} -> {chrom}:{pos} {ref}>{alt}")

    variant = genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases=alt,
    )

    # 预测区间
    sequence_length = 2048
    interval = variant.reference_interval.resize(sequence_length)

    # Scorer（你的原始设置）
    scorer = variant_scorers.CenterMaskScorer(
        width=None,
        aggregation_type=variant_scorers.AggregationType.DIFF_SUM_LOG2,
        requested_output=dna_client.OutputType.RNA_SEQ,
    )

    # 打分
    score_result = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=[scorer],
        organism=dna_client.Organism.HOMO_SAPIENS,
    )

    # 输出
    out_fn = f"{rsid}.txt"
    with open(out_fn, "w") as f:
        f.write(str(score_result[0].var))

    print(f"[OK] 保存结果到 {out_fn}")


def main():
    parser = argparse.ArgumentParser(description="Process all rsIDs from a VCF using AlphaGenome")
    parser.add_argument("--vcf", required=True, help="VCF file path")

    args = parser.parse_args()

    # 读取 rsIDs
    rsids = load_rsids_from_vcf(args.vcf)
    print(f"[INFO] 共找到 {len(rsids)} 个 rsID")

    # 初始化 AlphaGenome
    print("[INFO] 初始化 AlphaGenome 模型...")
    dna_model = dna_client.create(API_KEY)

    # 遍历每个 rsID —— 逐个处理
    for rsid in rsids:
        try:
            score_single_rsid(rsid, dna_model)
        except Exception as e:
            print(f"[ERROR] {rsid} 失败: {e}")


if __name__ == "__main__":
    main()
