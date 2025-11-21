import argparse
import gzip
import myvariant
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers

API_KEY = "AIzaSyD5Kht8QzCPkHeJ456_Tf_eBWirtKhmaRU"

def extract_rsids_from_csq(info_field):
    if "CSQ=" not in info_field:
        return None
    csq_data = info_field.split("CSQ=")[1]
    entries = csq_data.split(",")
    rsids = []
    for entry in entries:
        fields = entry.split("|")
        for f in fields:
            f = f.strip()
            if f.startswith("rs") and f[2:].isdigit():
                rsids.append(f)
    return rsids if rsids else None

def load_rsids_from_vcf(vcf_path):
    opener = gzip.open if vcf_path.endswith(".gz") else open
    all_rsids = []
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 8: continue
            info = parts[7]
            found = extract_rsids_from_csq(info)
            if found: all_rsids.extend(found)
    return list(set(all_rsids))

def rsid_to_variant_info(rsid, genome_build="hg19"):
    mv = myvariant.MyVariantInfo()
    out = mv.getvariant(rsid, fields="dbsnp")
    if not out or "dbsnp" not in out:
        raise ValueError(f"{rsid} 未找到 dbSNP 信息")
    dbsnp = out["dbsnp"]
    if dbsnp.get("vartype","").lower() != "snv":
        raise ValueError(f"{rsid} 不是 SNV")
    ref = dbsnp.get("ref")
    alt = dbsnp.get("alt")
    chrom = dbsnp.get("chrom")
    pos = dbsnp.get(genome_build, {}).get("start")
    if None in [ref, alt, chrom, pos]:
        raise ValueError(f"{rsid} 信息不完整")
    return f"chr{chrom}", pos, ref, alt

def score_single_variant(dna_model, rsid):
    try:
        chrom, pos, ref, alt = rsid_to_variant_info(rsid)
    except Exception as e:
        print(f"[WARN] 跳过 {rsid}: {e}")
        return None
    variant = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
    interval = variant.reference_interval.resize(2048)
    scorer = variant_scorers.CenterMaskScorer(
        width=None,
        aggregation_type=variant_scorers.AggregationType.DIFF_SUM_LOG2,
        requested_output=dna_client.OutputType.RNA_SEQ,
    )
    try:
        score_result = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=[scorer],
            organism=dna_client.Organism.HOMO_SAPIENS,
        )
        return chrom, pos, ref, alt, score_result[0].var
    except Exception as e:
        print(f"[ERROR] {rsid} 评分失败: {e}")
        return None

def main(vcf_path, output_csv="results.csv"):
    rsids = load_rsids_from_vcf(vcf_path)
    print(f"[INFO] 共找到 {len(rsids)} 个 rsID")
    dna_model = dna_client.create(API_KEY)
    results = []
    for rsid in rsids:
        print(f"[RUN] 处理 {rsid} ...")
        out = score_single_variant(dna_model, rsid)
        if out:
            chrom, pos, ref, alt, delta = out
            results.append({"rsid": rsid, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "delta_score": delta})
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"[FINISHED] 结果已写入 {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="批量解析 VCF 并打分")
    parser.add_argument("--vcf", type=str, required=True, help="VCF 文件路径")
    parser.add_argument("--out", type=str, default="results.csv", help="输出 CSV 文件")
    args = parser.parse_args()
    main(args.vcf, args.out)
