from matplotlib import pyplot as plt
def plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
                color_1, color_2, color_3, dot_style, label):
    plt.plot(epsA_calc, sigA_calc, color=color_1, linestyle=dot_style, label=label)
    plt.plot(epsR_calc, sigA_calc, color=color_2, linestyle=dot_style)
    plt.plot(epsV_calc, sigA_calc, color=color_3, linestyle=dot_style)
    plt.plot(epsA, sigA, color=color_1)
    plt.plot(epsR, sigA, color=color_2)
    plt.plot(epsV, sigA, color=color_3)
    #plt.show()
