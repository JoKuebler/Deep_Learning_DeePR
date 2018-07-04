class Result_Evaluator:

    @staticmethod
    def average_probability(final_props):

        if len(final_props) > 1:
            return sum(final_props) / len(final_props)
        else:
            return 0

    @staticmethod
    def false_negatives(results):

        counter = 0

        for result in results:
            if result < 0.6:
                counter += 1

        print('False Negatives', str(counter))
        return counter

    @staticmethod
    def false_positives(results):

        counter = 0

        for result in results:
            if result > 0.5:
                counter += 1

        print('False Positives', str(counter))
        return counter

    def sensitivity(self, results):

        print('Sensitivity')
        return (len(results) - self.false_negatives(results)) / len(results)

    def specificity(self, results):

        print('Specifity')
        return (len(results) - self.false_positives(results)) / len(results)

    def calculate_metrices(self, results, split_index):

        test_pos = results[0:split_index]
        test_neg = results[split_index + 1:]

        sensitivity = self.sensitivity(test_pos)
        specificity = self.specificity(test_neg)

        return [sensitivity, specificity]


