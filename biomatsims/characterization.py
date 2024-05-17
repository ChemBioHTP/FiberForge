



def calculate_var_over_time(xvg_file):
    with open(xvg_file, 'r') as f:
        lines = f.readlines()
        time_data = []
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                time_data.append(list(map(float, line.split())))
    return time_data