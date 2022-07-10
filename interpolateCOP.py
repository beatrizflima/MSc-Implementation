import csv
import math
import numpy as np
import pandas as pd

# Configuration
path_qh = r"pie_script\qh.csv"
path_qc = r"pie_script\qc.csv"
path_qh_needed = r"pie_script\qh_needed_reg.csv"
path_results = r"res.csv"
path_output_dataset = r"list_res.csv"
global_min = -1e100

# Code

def verify_input(qh, qc, qh_needed):
    if len(qh_needed) == 0:
        print("Empty qh_needed")
        exit(1)

    if len(qh) == 0 or len(qc) == 0:
        print("Empty dataset")
        exit(1)

    if len(qh) != len(qc):
        print("Missmatch in dimension [height]")
        exit(1)

    if len(qh[0]) != len(qc[0]):
        print("Missmatch in dimension [width]")
        exit(1)

# returns csv as a list of lists. Removes the header
def get_csv_rows(path):
    file = open(path)
    csv_reader = csv.reader(file)
    rows = []
    for row in csv_reader:
        row = [float(c) for c in row]
        rows.append(row)
    file.close()
    return rows

def store_list_of_list_as_csv(matrix, path):
    with open(path, 'w', newline='', encoding='utf-8') as file:
        csv_writer = csv.writer(file)
        for row in matrix:
            csv_writer.writerow(row)

def calc_qh_needed(qh_needed_path, qh_path):
    qh_needed = np.loadtxt(open(qh_needed_path,"rb"))
    '''qh_headers = np.loadtxt(open(qh_path,"rb"), delimiter=",")
    qh = qh_headers[:,1:len(qh_headers)]
    qh = qh[1:len(qh_headers),:]
    reg = math.ceil(np.max(qh_needed)/np.max(qh))

    for i in range(len(qh_needed)):
        if qh_needed[i] != 0.0:
            qh_needed[i] = qh_needed[i]/reg'''
    return qh_needed

def init_empty_matrix(rows, columns):
    matrix = []
    for row in range(rows):
        matrix.append([0] * columns)
    return matrix

def find_qh_interpolated_coords(qh, searched_val, nrows, ncolumns):
    coords_list = []

    # search by row
    for row in range(1, nrows):
        for col in range(1, ncolumns - 1):
            val1 = qh[row][col]
            val2 = qh[row][col + 1]
            if searched_val == val1:
                coords_list.append((col, row))
            elif (searched_val > val1 and searched_val < val2) or (searched_val < val1 and searched_val > val2):
                ratio = (searched_val - val1) / (val2 - val1)
                x_val = (1 - ratio) * (col) + ratio * (col + 1)
                coords_list.append((x_val, row))
        
        if searched_val == qh[row][-1]:
            coords_list.append((ncolumns - 1, row))
    
    # search by column
    for col in range(1, ncolumns):
        for row in range(1, nrows - 1):
            val1 = qh[row][col]
            val2 = qh[row + 1][col]
            if (searched_val > val1 and searched_val < val2) or (searched_val < val1 and searched_val > val2):
                ratio = (searched_val - val1) / (val2 - val1)
                y_val = (1 - ratio) * (row) + ratio * (row + 1)
                coords_list.append((col, y_val))
    
    return coords_list


def find_iterpolated_cop(cop, coords_list, nrows, ncolumns):
    max_cop = global_min
    freq_mf_pairs = []

    for coords in coords_list:
        x = coords[0]
        y = coords[1]
        x1 = math.floor(x)
        x2 = math.ceil(x)
        y1 = math.floor(y)
        y2 = math.ceil(y)

        if x1 == x2:
            x_ratio = 1
        else:
            x_ratio = (x - x1) / (x2 - x1)

        if y1 == y2:
            y_ratio = 1
        else:
            y_ratio = (y - y1) / (y2 - y1)
        
        point_cop1 = (1 - x_ratio) * cop[y1][x1] + x_ratio * cop[y2][x2]
        point_cop2 = (1 - y_ratio) * cop[y1][x1] + y_ratio * cop[y2][x2]
        mass_flow = (1 - x_ratio) * cop[0][x1] + x_ratio * cop[0][x2]
        freq = (1 - y_ratio) * cop[y1][0] + y_ratio * cop[y2][0]

        if isinstance(x, int):
            point_cop = point_cop2
        elif isinstance(y, int):
            point_cop = point_cop1
        #elif x==0:
            #point_cop = 0
        else:
            print("Both coords are not int")
            exit(1)

        
        if point_cop > max_cop:
            max_cop = point_cop
            freq_mf_pairs = [(freq, mass_flow)]
        elif point_cop == max_cop:
            freq_mf_pairs.append((freq, mass_flow))


    return max_cop, freq_mf_pairs


# Logic

def main():
    qh = get_csv_rows(path_qh)
    qc = get_csv_rows(path_qc)
    header = qh[0]
    qh_needed = calc_qh_needed(path_qh_needed, path_qh)
    #qh_needed = [14,0,15,0,16]
    verify_input(qh, qc, qh_needed)

    ncolumns = len(header)
    nrows = len(qh)
    result = init_empty_matrix(nrows, ncolumns)
    q_needed_and_coords = []

    for searched_val in qh_needed:
        q_needed_and_coords.append((searched_val, find_qh_interpolated_coords(qh, searched_val, nrows, ncolumns)))


    max_value = global_min
    idxs = []

    result[0] = qh[0] # first row is just copied over
    for row in range(1 ,nrows):
        result[row][0] = qh[row][0]     # first column is just copied over
        for col in range(1, ncolumns):
            qh_val = qh[row][col]
            qc_val = qc[row][col]
            result_cop = qh_val / (qh_val - qc_val)
            result[row][col] = result_cop

            if result_cop > max_value:
                max_value = result_cop
                idxs = [(row, col)]
            elif result_cop == max_value:
                idxs.append((row, col))
            
    store_list_of_list_as_csv(result, path_results)
    

    output_dataset = [["Q", "COP", "f", "mass"]]

    for qh_needed, coords_list in q_needed_and_coords:
        if qh_needed == 0:
            output_dataset.append([0, 0, 0, 0])
        #print(f"Qh needed: {qh_needed}")
        #print(f"\tMax COP: {round(max_cop, 5)}") 
        else:
            max_cop, freq_mf_pairs = find_iterpolated_cop(result, coords_list, nrows, ncolumns)
            max_cop = round(max_cop, 5)
            i = 0
            for freq_mf_pair in freq_mf_pairs:
                freq = round(freq_mf_pair[0], 5)
                mass_flow = round(freq_mf_pair[1], 5)
                #print(f"\t{i}) Freq: {freq}, Mass Flow: {mass_flow}")
                output_dataset.append([qh_needed, max_cop, freq, mass_flow])
                i += 1
        print()
        
    store_list_of_list_as_csv(output_dataset, path_output_dataset)
        
main()