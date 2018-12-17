import keras
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_score, recall_score, f1_score


class Histories(keras.callbacks.Callback):

    def __init__(self, out_path='/ebio/abt1_share/update_tprpred/code/src/network_data/', out_fn='best_model.h5'):
        self.f1 = 0
        self.path = out_path
        self.fn = out_fn
        self.threshold = 0.9

    def on_train_begin(self, logs={}):
        return

    def on_train_end(self, logs={}):
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):

        cv_pred = self.model.predict(self.validation_data[0])

        tp = cv_pred[:, 0] > self.threshold

        target = self.validation_data[1]

        one_dim_target = [0 if x[1] == 1.0 else 1 for x in target]

        f1 = 0

        tp_count = 0
        fp_count = 0
        fn_count = 0

        for index, elem in enumerate(tp):
            if tp[index] == one_dim_target[index] and elem:
                tp_count += 1
            elif tp[index] != one_dim_target[index] and not elem:
                fn_count += 1
            elif tp[index] != one_dim_target[index] and elem:
                fp_count += 1

        print('TP: ', tp_count)
        print('FN: ', fn_count)
        print('FP: ', fp_count)

        if tp_count > 0:
            f1 = (2 * tp_count) / (2 * tp_count + fn_count)
            prec_val = tp_count / (tp_count + fp_count)
            sens_val = tp_count / (tp_count + fn_count)

        if self.f1 < f1:
            f1_val2 = '%.3f' % f1
            sens_val2 = '%.3f' % sens_val
            prec_val2 = '%.3f' % prec_val
            print("Best F1 score: %s (prec: %s, sens: %s)" % (f1_val2, prec_val2, sens_val2))
            self.f1 = f1
            self.model.save(self.path + self.fn, overwrite=True)
            # Write model to json
            with open(self.path + 'best_model.json', 'w') as json_file:
                json_file.write(self.model.to_json())

        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return
