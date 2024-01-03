import numpy as np
def _find_closest_indices(lst, target_values):
        # Ensure lst is sorted for binary search
        sorted_indices = np.argsort(lst)
        sorted_lst = lst[sorted_indices]

        # Use binary search to find the closest indices
        closest_indices = np.searchsorted(sorted_lst, target_values, side='left')
        closest_indices = np.clip(closest_indices, 1, len(lst) - 1)  # Avoids index out of bounds

        # Determine the closer index based on absolute differences
        left_diff = np.abs(lst[sorted_indices[closest_indices - 1]] - target_values)
        right_diff = np.abs(lst[sorted_indices[closest_indices]] - target_values)

        closest_indices = np.where(left_diff < right_diff, closest_indices - 1, closest_indices)
        closest_indices = sorted_indices[closest_indices]

        return closest_indices

def test_find_closest_indices():
    # Sample data
    lst = np.array([1, 3, 5, 7, 9])

    # Target values
    target_values = np.array([1, 9, 5, 7])

    # Test the function
    closest_indices = _find_closest_indices(lst, target_values)

    print(closest_indices)
    # Expected results
    expected_indices = np.array([0, 4, 2, 3])

    # Check if the results match the expected values
    assert np.array_equal(closest_indices, expected_indices), "Test failed!"

test_find_closest_indices()