# filters out high p values for tprpred predictions
def filter_high_pvalues(input_file, output_directory):

    values = ['e+00', 'e-01', 'e-02', 'e-03']

    open_file = open(input_file, 'r')
    new_file = open(output_directory + 'pvalue_filtered' + '.txt', 'w')
    flag = False

    for line in open_file:

        pvalue_string = line[65:78]

        if flag:
            flag = False
            continue

        if any(x in pvalue_string for x in values):
            flag = True

        else:
            new_file.write(line)

    open_file.close()
    new_file.close()


if __name__ == '__main__':

    filter_high_pvalues('/ebio/abt1_share/update_tprpred/data/Pfam/tprpred_predictions.txt',
                        '/ebio/abt1_share/update_tprpred/data/Pfam/'
                        )
