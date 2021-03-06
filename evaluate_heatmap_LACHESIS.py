#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @File    : LA.py
# @Date    : 18-9-14
# @Author  : luyang(luyang@novogene.com)
import sys
import os
import argparse
from itertools import product
from math import ceil
import numpy as np
from scipy.stats import chi2_contingency
import xlsxwriter


def parse_input():
    parser = argparse.ArgumentParser(
        description='A Tool to Evaluate Heatmap generated by LACHESIS\nAuthour:luyang@novogene.com\n')
    parser.add_argument('--heatmap', '-m', required=True, dest="heatmap_file", default='heatmap.txt',
                        help='path to heatmap.txt')
    parser.add_argument('--break', '-b', required=True, dest="break_file", default='heatmap.chrom_breaks.txt',
                        help='path to heatmap.chrom_breaks.txt')
    parser.add_argument('--threshold', '-t', required=False, dest="p_value", default=0.0001, help='')
    parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    return args


def main():
    args = parse_input()
    path = os.getcwd()
    print(path)
    heatmap_file = args.heatmap_file
    break_file = args.break_file
    p_value = float(args.p_value)
    # if len(sys.argv) < 3:
    #     heatmap_file = 'heatmap.txt'
    #     break_file = 'heatmap.chrom_breaks.txt'
    #     p_value = 0.0001
    # else:
    #     heatmap_file = sys.argv[1]
    #     break_file = sys.argv[2]
    #     p_value = sys.argv[3]
    heatmap = {}
    chrom_breaks = []
    with open(break_file) as f:
        for line in f.readlines():
            chrom_breaks.append(ceil(float(line.strip())))
        tmp = []
        [tmp.append(i) for i in chrom_breaks if not i in tmp]
        chrom_breaks = tmp

    with open(heatmap_file) as f:
        next(f)
        for line in f.readlines():
            line = line.split()
            x = int(line[0])
            y = int(line[1])
            heat = float(line[2])
            if x in heatmap.keys():
                heatmap[x][y] = heat
            else:
                heatmap[x] = {}
                heatmap[x][y] = heat

    heatmap_array = np.zeros([len(heatmap), len(heatmap)], dtype=np.float64)
    for x in heatmap.keys():
        for y in heatmap[x].keys():
            heatmap_array[x, y] = heatmap[x][y]

        chrom_array = np.zeros([len(chrom_breaks) - 1, len(chrom_breaks) - 1])
    for i in range(len(chrom_breaks) - 1):
        for j in range(len(chrom_breaks) - 1):
            # x, dx = chrom_breaks[i:i + 2]
            # y, dy = chrom_breaks[j:j + 2]
            x, dx = chrom_breaks[i], chrom_breaks[i + 1] - chrom_breaks[i]
            y, dy = chrom_breaks[j], chrom_breaks[j + 1] - chrom_breaks[j]
            up_left = np.sum(heatmap_array[x:x + dx, y:y + dy])
            up_right = np.sum(heatmap_array[x:x + dx, :]) - up_left
            down_left = np.sum(heatmap_array[:, y:y + dy]) - up_left
            down_right = np.sum(heatmap_array) - up_left - up_right - down_left
            cross_table = np.array([
                [
                    up_left / dx / dy,
                    up_right / (dx) / (len(heatmap) - dy)
                ],
                [
                    down_left / (dy) / (len(heatmap) - dx),
                    down_right / (len(heatmap) * len(heatmap) - len(heatmap) * dy - len(heatmap) * dx + dx * dy)]
            ])
            chi2, p, dof, ex = chi2_contingency(cross_table)
            chrom_array[i, j] = p
    [rows, cols] = chrom_array.shape
    diag = np.diag(chrom_array)
    true_positive = 0
    false_negative = 0
    for i in range(len(chrom_array)):
        if diag[i] <= p_value:
            true_positive += 1
        else:
            false_negative += 1
    false_positive = 0
    true_negative = 0
    for i in range(rows):
        for j in range(cols):
            if chrom_array[i, j] <= p_value:
                false_positive += 1
            else:
                true_negative += 1
    false_positive = false_positive - true_positive
    total = np.array([
        [true_positive, false_positive],
        [false_negative, true_negative]
    ])
    recall = true_positive / (true_positive + false_negative)
    precision = true_positive / (true_positive + false_positive)
    f_score = 2 / (1 / recall + 1 / precision)
    # print('Recall:{:.2%}'.format(recall))
    # print('Precision:{:.2%}'.format(precision))
    print('F_score:{:.4f}'.format(f_score))
    outpath = path + '/heatmap.xlsx'
    workbook = xlsxwriter.Workbook(outpath)
    worksheet1 = workbook.add_worksheet('heatmap')
    worksheet2 = workbook.add_worksheet('score')
    scope = ''.join(['B2:', chr(ord('A') + rows), str(rows + 1)])
    total_fmt = workbook.add_format({'align': 'center', 'num_format': 0.00000})
    significant_fmt = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006', 'num_format': 0.00000})
    for i in range(rows + 1):
        for j in range(cols + 1):
            if i == 0 or j == 0:
                worksheet1.write_number(i, j, max(i, j))
            else:
                worksheet1.write_number(i, j, chrom_array[i - 1, j - 1], total_fmt)
    worksheet1.conditional_format(scope,
                                  {'type': 'cell', 'criteria': '<=', 'value': p_value, 'format': significant_fmt})
    worksheet2.write('D4', float('{:.4f}'.format(f_score)))
    for row, data in enumerate(total):
        worksheet2.write_row(row+1, 1, data)
    worksheet2.write_column(1,0,['Predicted condition positive', 'Predicted condition negative'])
    worksheet2.write_row(0,1,['Condition positive', 'Condition negative'])
    workbook.close()


if __name__ == "__main__":
    main()
