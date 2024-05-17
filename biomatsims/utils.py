import Bio
import Bio.SeqUtils

def single_to_triple(seq):
	seq = [Bio.SeqUtils.IUPACData.protein_letters_1to3[aa] for aa in seq]
	seq = ''.join(seq)
	t = iter(seq)
	trip_seq = ' '.join((a+b+c).upper() for a,b,c in zip(t,t,t))
	return trip_seq

aa_charge = {
    'A': 0,
    'C': 0,
    'D': -1,
    'E': -1,
    'F': 0,
    'G': 0,
    'H': 1,
    'I': 0,
    'K': 1,
    'L': 0,
    'M': 0,
    'N': 0,
    'P': 0,
    'Q': 0,
    'R': 1,
    'S': 0,
    'T': 0,
    'V': 0,
    'W': 0,
    'Y': 0
}


def search_pattern_in_consecutive_lines(file_path, pattern1, pattern2):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines) - 1):
            line1 = lines[i]
            line2 = lines[i + 1].strip()
            if pattern1 in line1 and pattern2 in line2:
                return i + 1
    
    return False

def add_line_at_line_number(file_path, line_number, line_content):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    lines.insert(line_number, line_content + '\n')

    with open(file_path, 'w') as file:
        file.writelines(lines)


if __name__ == "__main__":
	res = single_to_triple('GAG')
	print(res)
	 

