import numpy as np
import pandas as pd
from Bio import SeqIO
import sys

# --- 1. File Parsing Functions ---

def parse_fasta(file_path):
    """Parses a FASTA file and returns a list of headers (IDs) and sequences."""
    headers = []
    sequences = []
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            headers.append(record.id)
            sequences.append(str(record.seq).upper()) 
    except FileNotFoundError:
        print(f"Error: Input FASTA file not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing FASTA file: {e}", file=sys.stderr)
        sys.exit(1)
        
    if len(sequences) < 2:
        print(f"Error: Input file {file_path} must contain at least two sequences.", file=sys.stderr)
        sys.exit(1)

    return headers, sequences

def parse_scoring_matrix(file_path):
    """
    Parses a scoring matrix file (e.g., PAM, BLOSUM) into a dictionary of dictionaries.
    Uses robust pandas parsing to handle comments ('#') and mixed whitespace (sep='\s+').
    """
    try:
        # 這些檔案包含註解行 (#) 和不規則的空白
        df = pd.read_csv(
            file_path, 
            sep='\s+', # 處理一個或多個空白
            index_col=0, 
            engine='python', 
            comment='#' # 忽略註解行
        )
        
        df = df.dropna(axis=1, how='all')
        
        df.columns = [col.upper() for col in df.columns]
        df.index = [idx.upper() for idx in df.index]
        
        return df.to_dict()
    except FileNotFoundError:
        print(f"Error: Scoring matrix file not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing scoring matrix: {e}", file=sys.stderr)
        sys.exit(1)

# --- 2. Global Alignment (Needleman-Wunsch) ---
# 全域比對邏輯 (已通過測試)

def perform_global_alignment(seq1, seq2, score_matrix, gap):
    """
    Performs global alignment using the Needleman-Wunsch algorithm.
    """
    n = len(seq1)
    m = len(seq2)
    dp_table = np.zeros((n + 1, m + 1))
    trace_table = [[None for _ in range(m + 1)] for _ in range(n + 1)]

    # Initialization
    for i in range(1, n + 1):
        dp_table[i][0] = i * gap
        trace_table[i][0] = 'up'
    for j in range(1, m + 1):
        dp_table[0][j] = j * gap
        trace_table[0][j] = 'left'
    trace_table[0][0] = 'stop'

    # Filling
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            try:
                sub_score = score_matrix[seq1[i - 1]][seq2[j - 1]]
            except KeyError:
                sub_score = gap 
            
            match_score = dp_table[i - 1][j - 1] + sub_score
            delete_score = dp_table[i - 1][j] + gap
            insert_score = dp_table[i][j - 1] + gap

            max_score = max(match_score, delete_score, insert_score)
            dp_table[i][j] = max_score

            if max_score == match_score:
                trace_table[i][j] = 'diag'
            elif max_score == delete_score:
                trace_table[i][j] = 'up'
            else:
                trace_table[i][j] = 'left'

    # Traceback
    aln1, aln2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        move = trace_table[i][j]
        if move == 'diag':
            aln1 += seq1[i - 1]
            aln2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif move == 'up':
            aln1 += seq1[i - 1]
            aln2 += '-'
            i -= 1
        elif move == 'left':
            aln1 += '-'
            aln2 += seq2[j - 1]
            j -= 1
        elif move == 'stop':
            break
            
    return [(aln1[::-1], aln2[::-1])]

# --- 3. Local Alignment (Smith-Waterman) ---
# 區域比對邏輯 (包含錯誤修正)

def perform_local_alignment(seq1, seq2, score_matrix, gap):
    """
    Performs local alignment using Smith-Waterman. Finds all optimal alignments,
    filters by length, and sorts by string order (protein1 then protein2).
    """
    n = len(seq1)
    m = len(seq2)
    
    dp_table = np.zeros((n + 1, m + 1))
    trace_table = [[[] for _ in range(m + 1)] for _ in range(n + 1)] # 儲存多個路徑

    max_score = 0
    max_positions = [] 

    # Fill the DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            try:
                sub_score = score_matrix[seq1[i - 1]][seq2[j - 1]]
            except KeyError:
                sub_score = gap 
            
            match_score = dp_table[i - 1][j - 1] + sub_score
            delete_score = dp_table[i - 1][j] + gap
            insert_score = dp_table[i][j - 1] + gap
            
            current_score = max(0, match_score, delete_score, insert_score)
            dp_table[i][j] = current_score

            if current_score == 0:
                continue 
            
            # 儲存 *所有* 導致最高分的移動 (使用 if, if, if)
            if current_score == match_score:
                trace_table[i][j].append('diag')
            if current_score == delete_score:
                trace_table[i][j].append('up')
            if current_score == insert_score:
                trace_table[i][j].append('left')

            # Update max-scoring positions
            if current_score > max_score:
                max_score = current_score
                max_positions = [(i, j)]
            elif current_score == max_score:
                if (i, j) not in max_positions:
                    max_positions.append((i, j))

    if max_score == 0:
        return [] 

    # --- Traceback all optimal paths ---
    all_alignments = set() 
    
    def traceback_paths(i, j, current_aln1="", current_aln2=""):
        """Recursively traces back all paths from (i, j) to a score of 0."""
        
        # Base case: Reached a score of 0
        if dp_table[i][j] == 0:
            all_alignments.add((current_aln1[::-1], current_aln2[::-1]))
            return

        # --- 關鍵修正：使用 `if`, `if`, `if` 而不是 `if/elif` ---
        # 這確保如果 trace_table[i][j] 包含 ['diag', 'up']，
        # 程式會 *同時* 探索這兩條路徑。
        for move in trace_table[i][j]:
            if move == 'diag':
                traceback_paths(i - 1, j - 1, current_aln1 + seq1[i - 1], current_aln2 + seq2[j - 1])
            if move == 'up':
                traceback_paths(i - 1, j, current_aln1 + seq1[i - 1], current_aln2 + '-')
            if move == 'left':
                traceback_paths(i, j - 1, current_aln1 + '-', current_aln2 + seq2[j - 1])

    # Start traceback from ALL positions that hold the max_score
    for i, j in max_positions:
        traceback_paths(i, j)

    # --- Filter and Sort based on assignment rules ---
    
    results_list = list(all_alignments)
    if not results_list:
        return []
        
    # Rule 1: Filter by maximum length
    max_len = max(len(aln1) for aln1, aln2 in results_list)
    filtered_results = [(aln1, aln2) for aln1, aln2 in results_list if len(aln1) == max_len]
    
    # Rule 2: Sort by protein1 string, then protein2 string
    sorted_results = sorted(filtered_results, key=lambda x: (x[0], x[1]))
    
    return sorted_results

# --- 4. Output Writer ---

def write_fasta(output_path, headers, alignments):
    """
    Writes the list of alignments to an output file in FASTA format.
    """
    try:
        with open(output_path, 'w') as f:
            for aln1, aln2 in alignments:
                # 遵循輸出的範例格式
                f.write(f">{headers[0]}\n")
                f.write(f"{aln1}\n")
                f.write(f">{headers[1]}\n")
                f.write(f"{aln2}\n")
    except IOError as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

# --- 5. Main Function (Required API) ---

def alignment(input_path, score_path, output_path, aln, gap):
    """
    Main function to run the alignment process.
    """
    
    # 1. Parse Input Files
    headers, sequences = parse_fasta(input_path)
    seq1 = sequences[0]
    seq2 = sequences[1]
    output_headers = [headers[0], headers[1]] 
    
    score_matrix = parse_scoring_matrix(score_path)
    
    try:
        gap_penalty = int(gap)
        if gap_penalty > 0:
            gap_penalty = -gap_penalty
    except ValueError:
        print(f"Error: Gap score '{gap}' must be an integer.", file=sys.stderr)
        sys.exit(1)

    # 2. Perform Alignment
    final_alignments = []
    if aln.lower() == "global":
        final_alignments = perform_global_alignment(seq1, seq2, score_matrix, gap_penalty)
    elif aln.lower() == "local":
        final_alignments = perform_local_alignment(seq1, seq2, score_matrix, gap_penalty)
    else:
        print(f"Error: Unknown alignment type '{aln}'. Use 'global' or 'local'.", file=sys.stderr)
        sys.exit(1)
        
    # 3. Write Output
    write_fasta(output_path, output_headers, final_alignments)

# --- 6. Command-Line Execution ---

if __name__ == "__main__":
    
    # 處理旗標 (flag) 參數
    if len(sys.argv) == 11:
        try:
            arg_input = sys.argv[sys.argv.index('--input') + 1]
            arg_score = sys.argv[sys.argv.index('--score') + 1]
            arg_aln = sys.argv[sys.argv.index('-aln') + 1]
            arg_gap = sys.argv[sys.argv.index('--gap') + 1]
            arg_output = sys.argv[sys.argv.index('--output') + 1]
        except (ValueError, IndexError):
            print("Error: Could not parse all required command-line arguments. Check flag names and order.", file=sys.stderr)
            sys.exit(1)
    else:
        # 備用的位置 (positional) 參數處理
        if len(sys.argv) != 6:
             print("Usage: python hw3.py <input> <score> <output> <aln> <gap> OR use all flags.", file=sys.stderr)
             sys.exit(1)
        arg_input = sys.argv[1]
        arg_score = sys.argv[2]
        arg_output = sys.argv[3]
        arg_aln = sys.argv[4]
        arg_gap = sys.argv[5] 
        
    try:
        alignment(
            input_path=arg_input,
            score_path=arg_score,
            output_path=arg_output,
            aln=arg_aln,
            gap=arg_gap
        )
    except Exception as e:
        sys.exit(1)
