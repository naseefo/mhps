


# # Import the os module, for the os.walk function
# import os


# def find(name, path):
#     for root, dirs, files in os.walk(path):
#         if name in files:
#             return os.path.join(root, name)

# print(find('r_se21', 'data\earthquakes'))

from random import random, seed

def rand_no(x_lower, x_upper):
    '''
    Function to generate random number
    x_lower     :   lower bound value
    x_upper     :   upper bound value
    '''
    random_number = x_lower + (x_upper - x_lower)*random()
    return random_number

if __name__ == '__main__':
    
    seed_value = 70 # Seed value for random generator
    seed(seed_value)

    print('Printing a list of 10 random numbers\n')
    for i in range(10):
        sample_value = rand_no(30, 35)
        print("Number %d = %4.2f"%(i, sample_value))