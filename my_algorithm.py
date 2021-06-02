import pandas as pd
import numpy as np
import math
import sys


col_names = {"bumphunter": {"chr": "chr", "indexStart": "cstart", "indexEnd": "cend", "p.value": "p", "L": "len"},
             "dmrff": {"chr": "chr", "start": "cstart", "end": "cend", "p.value": "p", "n": "len"},
             "combined": {"chr": "chr", "start": "cstart", "end": "cend", "p": "p", "length": "len"}}

some_stats = []


def get_methyl_arr(chr_df, empty_val=0):
    """ Function makes needed list for comparing methylations, cytozines' values are constructed to be in one line
    :param chr_df: dataframe from specific chromosome
    :param empty_val: value for cytozine not belonging to region
    :return: dictionary where keys ar chromosomes and values are list with each cytozines' values
    """
    try:
        max_pos = max(list(chr_df["cend"]))
    except:
        max_pos = 10
    # print(max_pos, min(list(chr_df["cstart"])))
    chr_vals = [float(empty_val) for _ in range(int(max_pos) + 1)]  # that list containing all values
    lens = 0
    for i in range(chr_df.shape[0]):
        start = int(list(chr_df["cstart"])[i])
        end = int(list(chr_df["cend"])[i])
        val = float(list(chr_df["my_val"])[i])
        chr_vals[start:end] = [val for _ in range(end+1 - start)]
        lens += end - start
    some_stats.append(lens/max_pos)
    return chr_vals


def get_combo(bumps, empty_val=0, if_skip_chry=False):
    """ Function constructs combined dictionary from several different methylation dictionaries
    :param bumps: dataframes containg bumps
    :param empty_val: value for cytozine not belonging to region
    :param if_skip_chry: whether to not include chromosome y into calculations
    :return: similar structure dictionary, but with combined p values
    """
    # 1. iterating through each chromosome/key, take each of dictionaries' lists
    # 2. combine to one list (use 2nd method)
    # 3. make similar structure methylation dictionary
    print("Making combination")
    method_name = ""  # which method was used, bumphunter or dmrff
    # modifying dataframes to suit my needs
    for bmp in bumps:
        cols = bmp.columns
        for c in col_names:
            all_found = True
            for k in col_names[c]:
                if k != "p.value":
                    if k not in cols:
                        all_found = False
            if all_found:
                method_name = c
                break
        if method_name == "":
            raise Exception("Method not found")
        new_cols = []
        for k in bmp.columns:
            if k in col_names[method_name]:
                new_cols.append(col_names[method_name][k])
            else:
                new_cols.append(k)
        bmp.columns = new_cols
        # print(bmp.columns)
        # calculating values for every region
        sum_vals = sum([(-1) * math.log(float(pval), 10) for pval in bmp["p"]])
        vals = [v / sum_vals for v in [(-1) * math.log(float(pval), 10) for pval in bmp["p"]]]
        bmp["my_val"] = vals
        bmp.sort_values("chr", ascending=True)
    chroms = set(list(bumps[0]["chr"]))
    my_dict = {}  # that combined dictionary
    for c in chroms:  # technically all of them should have same keys (because its chromosomes)
        if c == "chrY":
            if if_skip_chry:
                continue
        print("Working with chromosome " + str(c))
        chr_lines = [get_methyl_arr(b.loc[b["chr"] == c], empty_val=empty_val) for b in bumps]
        max_len = max([len(ch) for ch in chr_lines])
        for ch in range(len(chr_lines)):
            for i in range(max_len - len(chr_lines[ch])):
                chr_lines[ch].append(empty_val)
        chr_vals = []  # that list that will have combined values
        for i in range(max_len):
            chr_vals.append(np.mean([chr_lines[j][i] for j in range(len(chr_lines))]))  # *(2**coef)
        my_dict[c] = chr_vals
    del chr_vals, bumps, chroms
    return my_dict


def compare_bumps(bumps, empty_val=0, if_print_every_chr=False, if_skip_chry=False):
    """ Make a dictionary that would be useful in comparing different bumps
    :param bumps: dataframe containg bumps
    :param empty_val: value for cytozine not belonging to region
    :param if_print_every_chr: whether to print result from every chromosome
    :param if_skip_chry: whether to not include chromosome y into calculations
    :return: similar structure dictionary but instead of p values it would be sth like confusion matrix
    """
    # Note: have different key for which file is which value
    # Note: have different keys for statistics (like overlap between each)
    # In short, make similar structure dictionary, but instead of p value use 2^n numbers where 2^0 for first file,
    # then 2^1 for second. What it means is that if specific cytozine belongs to region in fisrt and second, but not
    # third dictionary, it would have value 3, if to 2nd and 3rd, but not 1st, then 6 and so on
    print("Comparing two sets")
    # my_dict = {}  # that combined dictionary

    method_name = ""  # which method was used, bumphunter or dmrff
    for bmp in bumps:
        cols = bmp.columns
        for c in col_names:
            all_found = True
            for k in col_names[c]:
                if k != "p.value":
                    if k not in cols:
                        all_found = False
            if all_found:
                method_name = c
                break
        if method_name == "":
            raise Exception("Method not found")
        new_cols = []
        for k in bmp.columns:
            if k in col_names[method_name]:
                new_cols.append(col_names[method_name][k])
            else:
                new_cols.append(k)
        bmp.columns = new_cols
        vals = [1 for _ in range(len([1 for _ in range(len(bmp["p"]))]))]
        bmp["my_val"] = vals
        bmp.sort_values("chr", ascending=True)
    chroms = set(list(bumps[0]["chr"])) # technically all of them should have same keys (because its chromosomes)
    final_s = 0  # skaitiklis
    final_v = 0  # vardiklis
    for c in chroms:
        if c == "chrY":
            if if_skip_chry:
                continue
        print("Working with chromosome "+str(c))
        chr_lines = [get_methyl_arr(b.loc[b["chr"] == c], empty_val=empty_val) for b in bumps]
        max_len = max([len(ch) for ch in chr_lines])
        for ch in range(len(chr_lines)):
            for i in range(max_len - len(chr_lines[ch])):
                chr_lines[ch].append(empty_val)
        conf = [0, 0, 0, 0]  # 4 values: (both marked as regions), (1st marked but not 2nd), (vice versa), (none marked)
        v = 0  # vardiklis
        for i in range(max_len):
            vals = [chr_lines[j][i] for j in range(len(chr_lines))]
            if vals[0] > empty_val and vals[1] > empty_val:
                conf[0] += 1
                v += 1
            elif vals[0] > empty_val >= vals[1]:
                conf[1] += 1
                v += 1
            elif vals[0] <= empty_val < vals[1]:
                conf[2] += 1
                v += 1
            elif vals[0] <= empty_val and vals[1] <= empty_val:
                conf[3] += 1
        # my_dict[c] = conf
        if if_print_every_chr:
            print("Overlap for " + c + ": " + str(float(conf[0]) / v))
        final_s += conf[0]
        final_v += v
    if final_v == 0:
        raise Exception("Final v is 0. Consider decreasing empty_value.")
    print("###################################")
    print("Final overlap: " + str(float(final_s) / final_v))
    print("###################################")
    return None


def make_table(my_dict, cut_method="mean"):
    """ Makes dataframe from combined dictionary
    :param my_dict: dictionary that has chromosomes as keys and list of cytozines a values
    :param cut_method: what cutoff method to use, possible values: "mean" - uses average, "percX" - uses Xth percentile
    :return: dataframe with bumps
    """
    # 1. Take each chromosome one by one
    # 2. Look for regions using significance level
    # 3. Add region to table (construct equal length dictionary with column names)
    # 4. Convert dictionary (from 3, not the methyl one) to dataframe
    print("Making table of bumps")

    my_df = {"chr": [], "length": [], "start": [], "end": [], "p": []}
    all_vals = []
    for ch in my_dict:
        all_vals.extend(my_dict[ch])
    cut = None
    if cut_method == "mean":
        cut = np.mean(all_vals)
    if cut_method[:4] == "perc":
        # avg = np.mean(all_vals)
        percentile = int(cut_method[4:])
        all_vals = [v for v in all_vals if v > 0]
        # av_dict = pd.DataFrame({"vals": all_vals})
        cut = np.percentile(all_vals, percentile)
        # pd_vals = pd.Series(all_vals)
        # av_dict.to_csv("C:\\Users\\simki\\OneDrive\\Desktop\\Univero failai\\8 semestras\\Bakalaurinis\\all_vals.csv")
        # for i in range(100):
        #     if np.percentile(all_vals, i) <= avg <= np.percentile(all_vals, i+1):
        #         print("Vidurkis lygus procentilei "+str(i))
        # print(avg, cut)
    if cut is None:
        raise Exception(
            "Wrong cut method: " + str(cut_method) + ". Use one of these: \"mean\", \"percX\" (e.g. \"perc50\")")
    for ch in my_dict:
        arr = my_dict[ch]
        if_reg = False  # if its region
        start = 0  # index of start position
        # all_vals = []
        for i in range(len(arr)):
            v = arr[i]
            # all_vals.append(v)
            if v >= cut and not if_reg:
                if_reg = True
                start = i
            if v < cut and if_reg:
                if_reg = False
                my_df["chr"].append(ch)
                my_df["length"].append(int(i - start))
                my_df["start"].append(start)
                my_df["end"].append(int(i - 1))
                my_df["p"].append(0.01)  # just any value so I could compare
        if if_reg:  # if chromosome ended but region didn't
            my_df["chr"].append(ch)
            my_df["length"].append(int((i + 1) - start))
            my_df["start"].append(start)
            my_df["end"].append(int(i))
            my_df["p"].append(0.01)
    my_df = pd.DataFrame(my_df)
    return my_df


# Base plan:
# 1. Make this a script in a way that it would use command line arguments
# 2. Only files with bumps should be necessary
# 3. Other arguments might be: name_for_final_file, empty_val


files = []
# peaks = []
empty_val_par = 0
if_print_every_chr_par = 0  # should be either 1 or 0 and then would be made into boolean
cut_method_par = "mean"
what_to_do_par = 1  # default 1 - combine bumps, -c or 2 - compare bumps
skip_y_chrom_par = 0
savename = ""  # actually if default two file names would be concatenated
start_of_pars = False

args = sys.argv
for a in range(len(args)):
    if a == 0:
        continue
    arg = args[a]
    if arg[0] == "-":
        start_of_pars = True
        if arg == "-e":
            empty_val_par = float(args[a+1])
        if arg == "-p":
            try:
                if float(args[a+1]) > 0:
                    if_print_every_chr_par = True
            except Exception as e:
                if_print_every_chr_par = False
        if arg == "-m":
            cut_method_par = args[a+1]
        if arg == "-c":
            what_to_do_par = 2
        if arg == "-s":
            savename = args[a+1]
        if arg == "-y":
            skip_y_chrom_par = 1
    else:
        if not start_of_pars:
            files.append(arg)

if savename == "":
    for f in files:
        savename += (f.split("/|\\")[-1]).split(".")[0]
        savename += "_"
    savename += "comb.csv"

skip_y_chrom_par = True if skip_y_chrom_par == 1 else False

# for f in files:
#     file = pd.read_csv(f)
#     bump_dict = get_methyl_dict(file, empty_val=empty_val_par)
#     peaks.append(bump_dict)

if what_to_do_par == 1:
    # comb_dict = get_combo(peaks, empty_val=empty_val_par)
    comb_table = make_table(
        get_combo(
            [pd.read_csv(f) for f in files],
            empty_val=empty_val_par, if_skip_chry=skip_y_chrom_par), cut_method=cut_method_par)
    if savename[-4:] != ".csv":
        savename += ".csv"
    comb_table.to_csv(savename)
    del comb_table
if what_to_do_par == 2:
    compare_bumps(
        [pd.read_csv(f) for f in files],
        empty_val=empty_val_par, if_print_every_chr=if_print_every_chr_par, if_skip_chry=skip_y_chrom_par)
    # not sure if I even need such dict

# print(np.mean(some_stats))

del savename, files, some_stats

print("End of program")

# example run
# python my_algorithm.py bumps_98876_age.csv bumps_147430_age.csv -s test.csv
