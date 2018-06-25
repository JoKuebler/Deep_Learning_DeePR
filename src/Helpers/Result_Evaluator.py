class Result_Evaluator:

    @staticmethod
    def average_probability(final_props):

        if len(final_props) > 1:
            return sum(final_props) / len(final_props)
        else:
            return 0
