import math
import os
import re
import pandas as pd
import PySimpleGUI as sg
def calc_dh_ds(seq, c, na):
    # constants
    # constant [kcal/(K*mol)]
    a = -0.0108
    # gas constant [kcal/(K*mol)]
    r = 0.00199
    # absolute zero
    dh_dic = {
        'AA': -9.1, 'AT': -8.6, 'AG': -7.8, 'AC': -6.5,
        'TA': -6.0, 'TT': -9.1, 'TG': -5.8, 'TC': -5.6,
        'GA': -5.6, 'GT': -6.5, 'GG': -11.0, 'GC': -11.1,
        'CA': -5.8, 'CT': -7.8, 'CG': -11.9, 'CC': -11.0
    }
    ds_dic = {
        'AA': -0.0240, 'AT': -0.0239, 'AG': -0.0208, 'AC': -0.0173,
        'TA': -0.0169, 'TT': -0.0240, 'TG': -0.0129, 'TC': -0.0135,
        'GA': -0.0135, 'GT': -0.0173, 'GG': -0.0266, 'GC': -0.0267,
        'CA': -0.0129, 'CT': -0.0208, 'CG': -0.0278, 'CC': -0.0266
    }
    dh = 0
    ds = 0
    for i in range(len(seq)-1):
        pair = seq[i]+seq[i+1]
        if pair in dh_dic:
            dh += dh_dic[pair]
            ds += ds_dic[pair]
    tm = dh/(a+ds+r*math.log(c/4))- 273.15 + 16.6*math.log10(na)
    return [dh, ds, tm]

# GUI
sg.theme('DarkTeal7')
layout = [ [sg.Text('Select the target file.')],
            [sg.Text('file name', size=(15, 1)), sg.Input(), sg.FileBrowse('Select', key='inputFilePath')],
            [sg.Button('submit', key='submit'), sg.Button('exit', key='exit')],
            [sg.Output(size=(80,20))]]
window = sg.Window('Oligonucleotide Melting Temperature Calculator', layout)

while True:
    event, values = window.read()

    if event in [sg.WIN_CLOSED, 'exit']:
        break

    else:
        path = values['inputFilePath']
        if not re.match(r'.*\.csv', path):
            print('Please import CSV File.')
            continue
        if not os.path.exists(path):
            print('The file doesn\'t exitst.')
            continue
        data = pd.read_csv(path)
        l = len(data.index)
        nuc_list = ['A', 'T', 'G', 'C']
        # generate nuc pair
        seqs = data['sequence']
        cs = data['c(mol)']
        nas = data['na(mol/L)']
        dh_ds = pd.DataFrame(map(calc_dh_ds, seqs, cs, nas))
        data['dh'] = dh_ds[0]
        data['ds'] = dh_ds[1]
        data['tm'] = dh_ds[2]
        print(data)
        out_path = re.sub('.csv', '_out.csv', path)
        data.to_csv(out_path)
        print('The file is successfully saved: ' + out_path)

window.close
