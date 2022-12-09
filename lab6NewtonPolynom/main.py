import math
from math import pi, sin, cos
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from numpy import linspace

plt.style.use('seaborn-v0_8')
matplotlib.use('TkAgg')


def divided_diff(x, y):
    '''
    function to calculate the divided
    differences table
    '''
    n = len(y)
    coef = np.zeros([n, n])
    # the first column is y
    coef[:, 0] = y

    for j in range(1, n):
        for i in range(n - j):
            coef[i][j] = \
                (coef[i + 1][j - 1] - coef[i][j - 1]) / (x[i + j] - x[i])

    return coef


def newton_poly(coef, x_data, x):
    '''
    evaluate the newton polynomial
    at x
    '''
    n = len(x_data) - 1
    p = coef[n]
    for k in range(1, n + 1):
        p = coef[n - k] + (x - x_data[n - k]) * p
    return p


def calc_euclid_norm(ar1, ar2):
    error = 0
    for i in range(len(ar1)):
        error += (ar1[i] - ar2[i]) ** 2
    return math.sqrt(error)


def print_polynom_for_given_cnt_of_points(ddc):
    print("Newton`s polynom: ")
    print(ddc[0], end=' ')
    for i in range(1, len(ddc)):
        if ddc[i] >= 0:
            print('+', end=' ')
        else:
            print('-', end=' ')
        print(abs(ddc[i]), end='')
        for j in range(i):
            print('(x-x', j, ')', sep='', end='')
        print(end=' ')
    print()

    print("Newton`s polynom with given data: ")
    print(ddc[0], end=' ')
    for i in range(1, len(ddc)):
        if ddc[i] >= 0:
            print('+', end=' ')
        else:
            print('-', end=' ')
        print(abs(ddc[i]), end='')
        for j in range(i):
            if x_test[j] < 0:
                print('(x+', abs(x_test[j]), ')', sep='', end='')
            else:
                print('(x-', x_test[j], ')', sep='', end='')
        print(end=' ')


def do_things_for_given_functions():

    print_polynom_for_given_cnt_of_points(ddc)

    y_for_error = newton_poly(ddc, x_test, x)

    error = calc_euclid_norm(y, y_for_error)

    print("\nError:", math.sqrt(error))


print("Please input 1 for default task, 2 for y = -3x^3 - 5x^2 + 3x + 5, 3 for y = sin(x) + cos(x)^3")

type = int(input())

if type == 1:
    x_test = np.array([3, 4, 5, 6, 7])
    y_test = np.array([12, -1, 5, 8, 3])

    # get the divided difference coef
    ddc = divided_diff(x_test, y_test)[0, :]

    # evaluate on new data points
    x1 = linspace(2.5, 7.5, 1000)
    y_new = newton_poly(ddc, x_test, x1)
    y = y_test
    x = x_test
    cnt_points = 4
    do_things_for_given_functions()

    plt.figure(figsize=(10, 6))
    plt.xlabel('x')
    plt.ylabel('f(x)')

    plt.plot(x1, y_new, 'b')
    plt.plot(x_test, y_test, 'bo')

    plt.show()


elif type == 2:
    def fx1(x):
        return -3 * (x ** 3) - 5 * (x ** 2) + 3 * x + 5


    x_test0 = []
    print("Input number of test points:")
    cnt_points = int(input())
    cnt_points -= 1
    left = -2
    right = 2

    step = (right - left) / cnt_points

    while True:
        x_test0.append(left)
        left += step
        if left > right:
            break

    x_test = np.array(x_test0)
    y_test = np.array([fx1(x) for x in x_test])

    ddc = divided_diff(x_test, y_test)[0, :]

    x = linspace(-2, 2, 1000)
    y = -3 * (x ** 3) - 5 * (x ** 2) + 3 * x + 5
    y_new = newton_poly(ddc, x_test, x)

    do_things_for_given_functions()

    plt.xlabel = 'x'
    plt.ylabel = 'y = -3x^3 - 5x^2 + 3x + 5'
    plt.plot(x, y, 'r')
    plt.plot(x, y_new, 'b')
    plt.plot(x_test, y_test, 'bo')
    plt.show()

elif type == 3:
    def fx2(x):
        return sin(x) + cos(x) ** 3


    x = linspace(-pi, pi, 1000)
    y = [fx2(el) for el in x]

    x_test = []
    print("Input number of test points:")
    cnt_points = int(input())
    cnt_points -= 1
    left = -pi
    right = pi

    step = (right - left) / cnt_points

    while True:
        x_test.append(left)
        left += step
        if left > right:
            break

    x_test0 = np.array(x_test)

    x_test = x_test0
    y_test = np.array([fx2(x) for x in x_test])

    ddc = divided_diff(x_test, y_test)[0, :]
    y_new = newton_poly(ddc, x_test, x)

    do_things_for_given_functions()

    plt.xlabel = 'x'
    plt.ylabel = 'y = sin(x) + cox(x)^3'
    plt.plot(x, y, 'r')
    plt.plot(x, y_new, 'b')
    plt.plot(x_test, y_test, 'bo')
    plt.show()
