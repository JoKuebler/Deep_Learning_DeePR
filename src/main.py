import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def tprpred_plot(result_file):
    """
    Takes result file of TPRpred prediction and plots histogram
    :param result_file: TPRpred predictions including pvalues
    :return:
    """
    np_pval = []

    with open(result_file) as file:
        for line in file:
            p_val = line.split()[-1].split("=")[-1]
            np_pval.append(float(p_val))
            next(file)

    sns.distplot(np.array(np_pval))
    plt.show()


if __name__ == '__main__':

    tprpred_plot('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/tprpred_res.fa')



