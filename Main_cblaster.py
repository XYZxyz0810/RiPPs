import argparse
from Prodigal import run_prodigal # 从 Prodigal.py 中导入需要的函数
from Efetch_cblaster import run_fetch_nucl  # 从 Efetch2.py 导入 run_fetch_data 函数
from Hmmscan import run_hmmscan, process_hmmscan_output  # 假设 Hmmscan.py 中有 run_hmmscan_on_faa 函数
from Record_generate import process_prodigal_faa_output, comb_CDS_and_hmmscan_files, generate_cluster_prediction_tables, generate_precursor_fasta  # 假设 Record.py 中有 process_files 函数
from MMseqs import run_mmseqs2, parse_mmseqs_clusters, makelogo

def parse_arguments():
    """
    创建命令行解析器
    """
    parser = argparse.ArgumentParser(description="Run the bioinformatics pipeline.")
    parser.add_argument("-a", "--accession_file", help="Input file containing accession numbers")
    parser.add_argument("-u", "--updown_length", type=int, required=True, help="Upstream and downstream DNA length")
    parser.add_argument("pre_min", type=int, help="Minimal length of short peptide (integer).")
    parser.add_argument("pre_max", type=int, help="Maximal length of short peptide (integer).")
    parser.add_argument("out_dir", type=str, help="Output directory for results.")
    parser.add_argument("--nucl_file", type=str, help="Path to an existing nucleotide file.")
    parser.add_argument("--cblaster_file", type=str, help="Binary result file from cblaster.")
    return parser.parse_args()


def main():
    # 解析命令行参数
    args = parse_arguments()
    out_dir = args.out_dir  # 获取输出目录

    # 第一部分：获取核酸记录（使用 Efetch 获取数据）
    print("Running Nucl_Efetch...")
    if args.nucl_file:
        Nuclfile = args.nucl_file
        print(f"Using provided nucleotide file: {Nuclfile}")
    elif args.accession_file:
        print("Fetching data using accession file...")
        input_file = args.accession_file
        UDlength = args.updown_length
        Nuclfile = run_fetch_nucl(input_file, UDlength)
        print(f"Finished fetching data from accession. Data saved to {Nuclfile}")
    elif args.cblaster_file:
        print("Fetching data using cblaster binary result...")
        input_file = args.cblaster_file
        UDlength = args.updown_length
        Nuclfile = run_fetch_nucl(input_file, UDlength)
        print(f"Finished fetching data from cblaster result. Data saved to {Nuclfile}")
    else:
        raise ValueError("One of --nucl_file, --accession_file (-a), or --cblaster_file must be provided.")

    # 第二部分：处理核酸记录，生成faa文件和gbk文件（运行 Prodigal 预测 ORFs）
    print("Running Prodigal...")
    prodigal_faa, prodigal_gbk = run_prodigal(Nuclfile, out_dir)

    # 第三部分：处理faa输出文件，生成CDS记录文件和precursor记录文件
    print("Processing Prodigal output...")
    CDS_file, peptide_file = process_prodigal_faa_output(prodigal_faa, out_dir, args.pre_min, args.pre_max)
    print(f"Prodigal output processed.")

    # 第四部分：处理CDS记录文件，预测CDS的结构功能（运行 Hmmscan 并解析结果）
    print("Running Hmmscan on all predicted ORFs...")
    hmm_db_path = "Pfam-A.hmm"  # HMM 数据库路径
    hmmscan_output = run_hmmscan(prodigal_faa, hmm_db_path, out_dir)
    print("Parsing Hmmscan output...")
    parsed_hmmscan_results = process_hmmscan_output(hmmscan_output, out_dir)
    print(f"HMMScan parsing complete. Parsed results saved to {parsed_hmmscan_results}")

    # 第五部分：合并hmmscan结果文件和CDS预测结果文件，生成Cluster_prediction文件及可被MMseqs识别的多肽文件
    print("Processing final results...")
    # 处理文件并生成 CDS_prediction_file
    CDS_prediction_file = comb_CDS_and_hmmscan_files(out_dir, CDS_file, parsed_hmmscan_results)
    print(f"Results have been processed and saved to {CDS_prediction_file}")
    # 根据 accession 分组并生成最终分组表
    print("Generating grouped tables by accession...")
    Cluster_prediction_table = generate_cluster_prediction_tables(out_dir, CDS_prediction_file)
    print("Pipeline finished successfully!")
    # 处理文件peptide_file并生成被mmseqs2识别的peptide_output_file
    peptide_fasta_file = generate_precursor_fasta(out_dir, peptide_file)

    # 第六部分：处理多肽文件，使用mmseqs2聚类，并输出logo
    cluster_result_file = run_mmseqs2(peptide_fasta_file, out_dir)
    fasta_files = parse_mmseqs_clusters(cluster_result_file, peptide_fasta_file, out_dir, min_nodes=1)
    makelogo(fasta_files, out_dir)


if __name__ == "__main__":
    main()

#python3 Main_cblaster.py -a GT2_LanM_28.tsv -u 20000 20 150 GT2_LanM_28_out_dir
#python3 Main.py -a 1.txt -u 10000 20 150 out_dir --skip_fetch --nucl_file GT2_LanC_10kb_nucl.fasta