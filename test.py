# Your nested list
nested_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]

# Convert the nested list to a string
nested_list_str = '\n'.join('\t'.join(map(str, sublist)) for sublist in nested_list)
print(type(nested_list_str))

# Convert the string back to a nested list
restored_nested_list = [list(map(int, line.split('\t'))) for line in nested_list_str.split('\n')]
print(type(restored_nested_list))

# Print the original and restored nested lists
print("Original nested list:", nested_list)

print("Restored nested list:", restored_nested_list)

