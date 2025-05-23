with open('tvd.prop', 'w') as file:
    for i in range(1, 10000000):
        file.write(f"{i} 0\n")
