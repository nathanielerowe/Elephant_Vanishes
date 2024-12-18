class PROXMLMaker(object):
    def __init__(self):
        self.has_mode = False
        self.has_detector = False
        self.has_channel = False
        self.has_plotpot = False
        self._variations = []
        self.has_variations = False
        self.has_model = False

    def mode(self, m):
        self._mode = m
        self.has_mode = True
        return self

    def detector(self, d):
        self._detector = d
        self.has_detector = True
        return self

    def channel(self, name, unit, bins, true_min, true_max, true_nbin):
        self._channel_name = name
        self._channel_unit = unit
        self._channel_bins = bins
        self._channel_true_min = true_min
        self._channel_true_max = true_max
        self._channel_true_nbin = true_nbin
        self._subchannels = []
        self.has_channel = True
        return self

    def add_subchannel(self, name, plotname, color):
        assert(self.has_channel)
        self._subchannels.append(dict(name=name, plotname=plotname, color=color))
        return self

    def plotpot(self, p):
        self._plotpot = p
        self.has_plotpot = True
        return self

    def add_variation(self, typ, name, allow=True):
        if allow:
             self._variations["allowlist"] = dict(type=typ, name=name)
        else:
             self._variations["denylist"] = dict(type=typ, name=name)
        self.has_variations = True
        return self

    def model(self, tag):
        self._model_tag = tag
        self.has_model = True
        self._model_rules = []
        return self

    def add_model_rule(self, name, index):
        assert(self.has_model)
        self._model_rules.append(dict(index=index, name=name))
        return self
