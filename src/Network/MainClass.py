from src.Helpers.Result_Evaluator import Result_Evaluator
from src.Network.FeatureSelection import FeatureSelection
from src.Network.PredictorMain import MainPredictor

if __name__ == '__main__':

    # ------------ Initializing Objects -------------
    # Initalize Main Runner Object
    main_runner = MainPredictor()
    # Initialize necessary objects
    # Preprocessor with input data to predict
    Preprocessor = main_runner.init_preprocessor(main_runner.test_two)
    # Network
    Network = main_runner.init_network()
    # Feature Selector
    FeatureSelector = FeatureSelection()
    # Result Evaluator
    Result_Evaluator = Result_Evaluator()
    # -----------------------------------------------

    # --------------- Training ----------------------
    # Initialize objects when training model
    Training_data = main_runner.init_training_data(Preprocessor)
    # Traing network on data and apply cross validation
    main_runner.train_network(Network, Training_data)
    # main_runner.cross_validate(Network, Training_data)
    # -----------------------------------------------

    # --------------- Safe & Load -------------------
    # If model has to be saved after training
    # main_runner.save_network(Network)

    # Load network from saved json
    # main_runner.load_network(Network)
    # -----------------------------------------------

    # --------------- Prediction --------------------
    # Predict raw training data (fragments)
    # main_runner.predict_training_data(Network, Preprocessor_two, 0.8)

    # predict desired data with threshold (sequences)
    # main_runner.single_predict(Network, Preprocessor, 0.93)
    # -----------------------------------------------

    # --------------Evaluation ----------------------
    # Select Features
    # FeatureSelector.fit(Training_data[0], Training_data[1])

    # Predict Test Data and calculate Sensitivity and Specifity
    # results = Network.predict(Training_data[2])
    # print(Result_Evaluator.calculate_metrices(results, 7757))
    # ------------------------------------------------


























    # ref_support = [False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, True, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, True, True, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False, False, False, False, False,
    #                False, False, False, False, False, False, False, False, False, False]
    #
    # ref_ranking = [669, 299, 343, 434, 670, 567, 782, 212, 772, 260, 708, 184, 70, 578, 433, 579, 39, 836, 693, 765,
    #                707, 771, 295, 781, 102, 449, 668, 548, 681, 680, 637, 756, 221, 436, 198, 157, 809, 522, 682, 241,
    #                525, 524, 523, 683, 186, 298, 448, 547, 240, 164, 316, 605, 319, 811, 320, 142, 130, 806, 317, 239,
    #                424, 763, 57, 847, 422, 839, 315, 792, 562, 602, 318, 778, 526, 423, 69, 96, 44, 497, 519, 470, 606,
    #                329, 268, 97, 323, 128, 88, 271, 199, 269, 45, 734, 521, 43, 38, 270, 328, 498, 520, 95, 770, 242,
    #                375, 769, 594, 285, 776, 193, 753, 722, 131, 123, 367, 172, 217, 663, 661, 398, 256, 549, 619, 662,
    #                634, 775, 163, 626, 629, 843, 360, 381, 720, 303, 48, 744, 442, 574, 15, 301, 156, 208, 627, 330,
    #                290, 408, 628, 300, 704, 705, 302, 799, 63, 247, 110, 504, 288, 245, 503, 506, 192, 553, 505, 190,
    #                152, 151, 532, 24, 531, 150, 655, 149, 58, 552, 246, 630, 191, 546, 455, 544, 545, 543, 53, 667, 331,
    #                399, 542, 6, 144, 7, 666, 333, 534, 797, 195, 814, 813, 698, 178, 145, 332, 177, 321, 515, 804, 464,
    #                844, 430, 826, 426, 825, 824, 79, 21, 259, 363, 822, 823, 407, 362, 427, 267, 428, 692, 516, 258,
    #                429, 702, 253, 643, 644, 817, 162, 425, 701, 250, 160, 40, 337, 338, 768, 249, 528, 209, 336, 175,
    #                631, 161, 642, 527, 42, 41, 342, 111, 686, 89, 601, 365, 159, 687, 374, 599, 33, 36, 640, 35, 37, 80,
    #                373, 372, 16, 833, 364, 807, 600, 685, 664, 379, 726, 632, 789, 724, 550, 229, 289, 633, 621, 62, 25,
    #                225, 189, 397, 837, 378, 230, 226, 723, 227, 725, 620, 228, 757, 314, 709, 116, 388, 741, 384, 140,
    #                803, 382, 277, 108, 284, 810, 282, 386, 529, 735, 385, 593, 387, 800, 665, 530, 49, 383, 234, 252,
    #                101, 608, 609, 170, 418, 654, 218, 129, 291, 456, 61, 60, 752, 233, 417, 352, 294, 648, 52, 733, 419,
    #                51, 656, 340, 410, 147, 646, 784, 645, 148, 414, 412, 514, 9, 411, 509, 511, 512, 510, 415, 508, 513,
    #                76, 409, 413, 658, 146, 8, 176, 276, 275, 181, 576, 575, 82, 183, 173, 729, 50, 450, 83, 122, 180,
    #                623, 451, 179, 598, 622, 684, 577, 182, 612, 3, 831, 185, 802, 596, 801, 840, 447, 710, 774, 107,
    #                104, 435, 278, 103, 84, 796, 279, 377, 238, 350, 106, 773, 105, 46, 446, 313, 795, 794, 476, 580,
    #                474, 477, 139, 478, 99, 138, 136, 793, 653, 828, 563, 480, 475, 222, 187, 762, 479, 473, 137, 5, 746,
    #                248, 232, 758, 551, 392, 394, 495, 393, 280, 700, 570, 100, 283, 569, 496, 127, 126, 568, 73, 821, 4,
    #                75, 74, 494, 589, 590, 573, 588, 591, 64, 754, 134, 592, 188, 755, 788, 1, 273, 572, 571, 777, 322,
    #                736, 380, 13, 14, 12, 28, 272, 396, 261, 714, 203, 395, 366, 783, 204, 715, 507, 202, 716, 255, 165,
    #                81, 155, 566, 731, 717, 55, 370, 740, 17, 205, 730, 728, 829, 561, 469, 558, 90, 201, 557, 251, 556,
    #                554, 91, 732, 635, 560, 292, 293, 223, 466, 78, 727, 555, 92, 559, 30, 153, 174, 77, 310, 737, 309,
    #                659, 441, 660, 518, 56, 206, 453, 306, 308, 816, 324, 815, 171, 618, 307, 452, 305, 616, 617, 254,
    #                98, 281, 471, 500, 502, 143, 764, 348, 467, 216, 215, 818, 11, 18, 349, 200, 347, 211, 706, 501, 10,
    #                641, 499, 66, 311, 615, 458, 748, 420, 213, 244, 747, 312, 416, 359, 406, 32, 339, 325, 421, 834,
    #                636, 827, 224, 639, 457, 31, 20, 124, 118, 432, 54, 326, 846, 718, 838, 327, 132, 719, 26, 405, 266,
    #                264, 812, 657, 779, 265, 805, 790, 472, 431, 59, 27, 133, 760, 465, 460, 461, 743, 114, 647, 459,
    #                197, 462, 820, 220, 34, 296, 835, 463, 517, 117, 759, 219, 610, 742, 611, 1, 1, 690, 787, 688, 689,
    #                786, 65, 672, 355, 445, 443, 235, 356, 785, 358, 357, 109, 691, 135, 671, 47, 2, 214, 613, 444, 614,
    #                491, 115, 439, 112, 492, 353, 490, 113, 493, 376, 120, 780, 597, 624, 625, 438, 437, 125, 210, 391,
    #                703, 440, 354, 119, 19, 196, 538, 536, 341, 335, 652, 751, 830, 535, 389, 287, 344, 297, 361, 68, 67,
    #                750, 345, 848, 346, 262, 651, 650, 263, 537, 713, 304, 486, 141, 483, 489, 286, 841, 484, 154, 71,
    #                485, 712, 274, 482, 481, 231, 93, 94, 488, 711, 607, 749, 72, 487, 402, 808, 583, 158, 584, 581, 582,
    #                166, 403, 638, 167, 23, 390, 404, 738, 739, 400, 351, 169, 721, 168, 401, 595, 798, 29, 371, 679,
    #                454, 678, 586, 675, 673, 676, 649, 368, 236, 745, 565, 674, 564, 819, 585, 257, 369, 832, 237, 587,
    #                677, 194, 22, 533, 121, 604, 85, 541, 86, 694, 87, 761, 207, 243, 539, 468, 334, 767, 845, 699, 842,
    #                697, 766, 791, 603, 695, 696, 540]
    #
    # print(len(ref_ranking))
    # print(len(ref_support))
    # index = 0
    # while index < len(ref_ranking):
    #     print(ref_ranking[index:index+25])
    #     index += 25