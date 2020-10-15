from utils import general_pr_box


def test_general_pr_box_two_input_two_output_case():
    # x = y = 0 --> symmetric
    assert general_pr_box(0, 0, 0, 0) == 1 / 2
    assert general_pr_box(0, 1, 0, 0) == 0
    assert general_pr_box(1, 0, 0, 0) == 0
    assert general_pr_box(1, 1, 0, 0) == 1 / 2
    # x = 0 and y = 1 --> symmetric
    assert general_pr_box(0, 0, 0, 1) == 1 / 2
    assert general_pr_box(0, 1, 0, 1) == 0
    assert general_pr_box(1, 0, 0, 1) == 0
    assert general_pr_box(1, 1, 0, 1) == 1 / 2
    # x = 1 and y = 0 --> symmetric
    assert general_pr_box(0, 0, 1, 0) == 1 / 2
    assert general_pr_box(0, 1, 1, 0) == 0
    assert general_pr_box(1, 0, 1, 0) == 0
    assert general_pr_box(1, 1, 1, 0) == 1 / 2
    # x = 1 and y = 1 --> anti-symmetric
    assert general_pr_box(0, 0, 1, 1) == 0
    assert general_pr_box(0, 1, 1, 1) == 1 / 2
    assert general_pr_box(1, 0, 1, 1) == 1 / 2
    assert general_pr_box(1, 1, 1, 1) == 0
    return True


if __name__ == '__main__':
    test_general_pr_box_two_input_two_output_case()

