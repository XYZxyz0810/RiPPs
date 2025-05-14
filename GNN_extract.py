import pandas as pd
import re
from datetime import datetime


def clean_and_tail_gnn(input_file='GNN.csv', output_file='Temp_GNN.csv'):
    """
    Clean and preprocess the input GNN data file.

    Args:
        input_file (str): Path to the input CSV file. Defaults to 'GNN.csv'.
        output_file (str): Path to save the cleaned output file. Defaults to 'Temp_GNN.csv'.
    """
    # Read input file with GBK encoding
    df = pd.read_csv(input_file, encoding="gbk")

    # Remove rows with missing values in key columns
    df = df.dropna(subset=['SSN Cluster Number', 'Co-occurrence', 'Co-occurrence Ratio', 'Pfam'])

    # Convert SSN Cluster Number to integer, fill NA with 0
    df['SSN Cluster Number'] = pd.to_numeric(df['SSN Cluster Number'], errors='coerce').fillna(0).astype(int)

    # Convert Co-occurrence to float and round to 2 decimal places
    df['Co-occurrence'] = pd.to_numeric(df['Co-occurrence'], errors='coerce').round(2)

    def clean_ratio(val):
        """
        Clean and standardize the Co-occurrence Ratio values.

        Args:
            val: Input value to be cleaned

        Returns:
            str: Cleaned value in MM/DD format or original string if no date found
        """
        if pd.isna(val):
            return ""
        val = str(val)

        # Handle date formats: YYYY-MM-DD or YYYY/MM/DD
        date_match = re.search(r'(\d{4})[-/](\d{1,2})[-/](\d{1,2})', val)
        if date_match:
            month = int(date_match.group(2))
            day = int(date_match.group(3))
            return f"{month}/{day}"

        # Handle Chinese date format (e.g., "1月6日")
        cn_date_match = re.search(r'(\d{1,2})月(\d{1,2})日', val)
        if cn_date_match:
            month = int(cn_date_match.group(1))
            day = int(cn_date_match.group(2))
            return f"{month}/{day}"

        return val

    # Apply cleaning function to Co-occurrence Ratio column
    df['Co-occurrence Ratio'] = df['Co-occurrence Ratio'].apply(clean_ratio).astype(str)

    # Save cleaned data to output file
    df.to_csv(output_file, index=False, encoding='utf-8')
    print(f"[Step 1] Cleaned data saved to {output_file}")


def simplify_temp_file(input_file='Temp_GNN.csv', output_file='Temp_GNN_simplified.csv'):
    """
    Create a simplified version of the cleaned data with only essential columns.

    Args:
        input_file (str): Path to the cleaned input file. Defaults to 'Temp_GNN.csv'.
        output_file (str): Path to save the simplified output. Defaults to 'Temp_GNN_simplified.csv'.
    """
    df = pd.read_csv(input_file)

    # Select only relevant columns for further processing
    essential_columns = [
        'SSN Cluster Number',
        'Co-occurrence',
        'Co-occurrence Ratio',
        'Pfam',
        'Pfam Description'
    ]
    df_simplified = df[essential_columns]

    # Save simplified data
    df_simplified.to_csv(output_file, index=False)
    print(f"[Step 2] Simplified data saved to {output_file}")


def extract_pfam_matches(data_file='Temp_GNN_simplified.csv',
                         pfam_file='pfam_list.txt',
                         output_all='GNN_extract.txt',
                         output_filtered='GNN_extract_filtered.txt'):
    """
    Extract and aggregate Pfam matches from the simplified data.

    Args:
        data_file (str): Path to simplified data file. Defaults to 'Temp_GNN_simplified.csv'.
        pfam_file (str): Path to file containing Pfam IDs to match. Defaults to 'pfam_list.txt'.
        output_all (str): Path to save all matches. Defaults to 'GNN_extract.txt'.
        output_filtered (str): Path to save filtered matches (Co-occurrence > 0.2). Defaults to 'GNN_extract_filtered.txt'.
    """
    import collections

    # Load data and Pfam list
    df = pd.read_csv(data_file)
    with open(pfam_file, 'r') as f:
        pfam_list = [line.strip() for line in f if line.strip()]

    with open(output_all, 'w') as out_all, open(output_filtered, 'w') as out_filtered:
        for pfam in pfam_list:
            # Find rows containing the current Pfam ID
            matches = df[df['Pfam'].astype(str).str.contains(pfam, na=False)].copy()
            if matches.empty:
                continue

            # Initialize dictionary to aggregate ratio data by cluster
            ratio_data = collections.defaultdict(lambda: {
                'cooccur_sum': 0.0,
                'num_sum': 0,
                'den': None
            })

            # Process each matching row
            for _, row in matches.iterrows():
                cluster = int(row['SSN Cluster Number'])
                cooccur = float(row['Co-occurrence'])

                # Parse numerator and denominator from Co-occurrence Ratio
                ratio_str = str(row['Co-occurrence Ratio'])
                ratio_match = re.match(r'(\d+)/(\d+)', ratio_str)
                if ratio_match:
                    num = int(ratio_match.group(1))
                    den = int(ratio_match.group(2))
                else:
                    num, den = 0, None  # Default values if parsing fails

                # Initialize or validate denominator consistency
                if ratio_data[cluster]['den'] is None:
                    ratio_data[cluster]['den'] = den
                elif den != ratio_data[cluster]['den']:
                    print(
                        f"[Warning] Inconsistent denominator for cluster #{cluster} in Pfam {pfam}, skipping ratio aggregation.")
                    continue

                # Aggregate values
                ratio_data[cluster]['cooccur_sum'] += cooccur
                ratio_data[cluster]['num_sum'] += num

            # Sort clusters by number
            sorted_clusters = sorted(ratio_data.items())

            # 写入所有匹配结果（每个 Pfam 两行 + 空行）
            out_all.write(f"{pfam}\n")
            out_all.write(', '.join(
                f"#{cluster}({data['num_sum']}/{data['den']})" if data['den'] is not None
                else f"#{cluster}(NA)"
                for cluster, data in sorted_clusters
            ) + '\n')
            out_all.write('\n')  # 添加空行

            # 写入过滤后的结果（仅 cooccur > 0.2 且有 den）
            filtered_matches = ', '.join(
                f"#{cluster}({data['num_sum']}/{data['den']})"
                for cluster, data in sorted_clusters
                if data['cooccur_sum'] > 0.2 and data['den'] is not None
            )

            if filtered_matches:
                out_filtered.write(f"{pfam}\n")
                out_filtered.write(filtered_matches + '\n')
                out_filtered.write('\n')  # 添加空行

    print(
        f"[Step 3] Extraction complete. Results saved to: {output_all} (all matches) and {output_filtered} (filtered matches)")


if __name__ == "__main__":
    # Execute processing pipeline
    clean_and_tail_gnn()  # Step 1: Data cleaning
    simplify_temp_file()  # Step 2: Data simplification
    extract_pfam_matches()  # Step 3: Pfam extraction and analysis