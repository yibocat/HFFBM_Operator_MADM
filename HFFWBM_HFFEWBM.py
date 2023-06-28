import mohusets.fuzzynumbers as fns


def HFFWBM(alter, qrung, weight=None, p=1, q=1):
    def bonferroni_element(d1, d2, p_, q_):
        s1 = fns.algeb_multiply(d1.algeb_power(p_), d2.algeb_power(q_))
        s2 = fns.algeb_multiply(d1.algeb_power(q_), d2.algeb_power(p_))
        return fns.algeb_plus(s1, s2)

    dlist = []
    if weight is None:
        dl = alter
        for i, d1 in enumerate(dl):
            for j, d2 in enumerate(dl[:i]):
                dlist.append(bonferroni_element(d1, d2, p, q))
    else:
        dl = []
        assert len(alter) == len(weight), 'The number of alter is not equal to the number of weight!'
        for i in range(len(alter)):
            dl.append(alter[i].algeb_times(weight[i]))
        for i, d1 in enumerate(dl):
            for j, d2 in enumerate(dl[:i]):
                dlist.append(bonferroni_element(d1, d2, p, q))

    agge = fns.qrungdhfe(qrung, [0], [1])
    for agg in dlist:
        agge = fns.algeb_plus(agg, agge)
    agge = (agge.algeb_times(1 / (len(alter) * (len(alter) - 1)))).algeb_power(1 / (p + q))
    return agge


def HFFEWBM(alter, qrung, weight=None, p=1, q=1):
    def bonferroni_element(d1, d2, p_, q_):
        s1 = fns.eins_multiply(d1.eins_power(p_), d2.eins_power(q_))
        s2 = fns.eins_multiply(d1.eins_power(q_), d2.eins_power(p_))
        return fns.eins_plus(s1, s2)

    dlist = []
    if weight is None:
        dl = alter
        for i, d1 in enumerate(dl):
            for j, d2 in enumerate(dl[:i]):
                dlist.append(bonferroni_element(d1, d2, p, q))
    else:
        dl = []
        assert len(alter) == len(weight), 'The number of alter is not equal to the number of weight!'
        for i in range(len(alter)):
            dl.append(alter[i].eins_times(weight[i]))
        for i, d1 in enumerate(dl):
            for j, d2 in enumerate(dl[:i]):
                dlist.append(bonferroni_element(d1, d2, p, q))

    agge = fns.qrungdhfe(qrung, [0], [1])
    for agg in dlist:
        agge = fns.eins_plus(agg, agge)
    agge = (agge.eins_times(1 / (len(alter) * (len(alter) - 1)))).eins_power(1 / (p + q))
    return agge
