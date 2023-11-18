if __name__ == '__main__':
    l_fr0 = 0.1
    hx = 0.005
    lag = 0.01
    Nfr_with_fluid = int((l_fr0 - lag) / hx)+1
    Pinj = 100*10**5
    Pfr = [Pinj / 2 for i in range(Nfr_with_fluid)]
def calcK1c(l_fr0, hx, Nfr_with_fluid, Pfr):
    x_list = [hx * i for i in range(Nfr_with_fluid)]
    Sh = 5*10**5
    func_list = [(Pfr[i]-Sh) * ((x_list[i]+l_fr0)/(l_fr0-x_list[i])) ** 0.5 for i in range(Nfr_with_fluid)]
    integral_func = [(func_list[i]+func_list[i+1])/2*hx for i in range(Nfr_with_fluid-1)]
    K1 = 2/(3.14*l_fr0)**0.5*sum(integral_func)
    print(Nfr_with_fluid)
    print(x_list)
    print(K1)
    return K1