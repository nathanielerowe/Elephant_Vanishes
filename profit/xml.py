import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

@dataclass
class PROXMLSubChannel:
   name: str
   plotname: str = ""
   color: str = ""

   def to_xml(self, element):
       s = ET.SubElement(element, "subchannel")
       s.set("name", self.name)
       s.set("plotname", self.plotname)
       s.set("color", self.color)

@dataclass
class PROXMLChannel:
    name: str
    unit: str = ""
    bins: list[float] = field(default_factory=list)
    truebin_min: float = 0
    truebin_max: float = 2
    truebin_nbin: int = 200
    subchannel: list[PROXMLSubChannel] = field(default_factory=list)

    def to_xml(self):
        elem = ET.Element("channel")
        elem.set("name", self.name)
        elem.set("unit", self.unit)

        bins = ET.SubElement(elem, "bins")
        bins.set("edges", " ".join([str(b) for b in self.bins]))

        truebins = ET.SubElement(elem, "truebins")
        truebins.set("min", str(self.truebin_min))
        truebins.set("max", str(self.truebin_max))
        truebins.set("nbins", str(self.truebin_nbin))

        for s in self.subchannel:
            s.to_xml(elem)

        return elem

@dataclass
class PROXMLVariation:
    name: str
    type: str

    def to_xml(self, element, vartype):
        s = ET.SubElement(element, "%slist" % vartype)
        s.set("type", self.type)
        s.text = self.name

@dataclass
class PROXMLDetector:
    name: str

    def to_xml(self):
        e = ET.Element("detector")
        e.set("name", self.name)
        return e

@dataclass
class PROXMLModelRule:
    name: str
    index: int

    def to_xml(self, element):
        s = ET.SubElement(element, "rule")
        s.set("name", self.name)
        s.set("index", str(self.index))

@dataclass
class PROXMLModel:
    tag: str
    rule: list[PROXMLModelRule]

    def to_xml(self):
        e = ET.Element("model")
        e.set("tag", self.tag)
        for r in self.rule:
            r.to_xml(e)
        return e

@dataclass
class PROXMLFriendTree:
    filename: str
    treename: str

    def to_xml(self, element):
        s = ET.SubElement(element, "friend")
        s.set("filename", self.filename)
        s.set("treename", self.treename)

@dataclass
class PROXMLBranch:
    associated_subchannel: str
    name: str
    true_param_name: str
    true_L_name: str
    pdg_name: str
    additional_weight: str
    model_rule: str = "1"

    @staticmethod
    def get_associated_subchannel(mode, detector, channel, subchannel):
        return mode + "_" + detector + "_" + channel.name + "_" + subchannel.name

    def to_xml(self, element):
        s = ET.SubElement(element, "branch")
        s.set("associated_subchannel", self.associated_subchannel)
        s.set("name", self.name)
        s.set("true_param_name", self.true_param_name)
        s.set("true_L_name", self.true_L_name)
        s.set("pdg_name", self.pdg_name)
        s.set("additional_weight", self.additional_weight)
        s.set("model_rule", self.model_rule)

@dataclass
class PROXMLMCFile:
    filename: str
    treename: str
    pot: float
    branch: list[PROXMLBranch] = field(default_factory=list)
    scale: float = 1.
    maxevents: int = -1
    friend: list[PROXMLFriendTree] = field(default_factory=list)

    def to_xml(self):
        e = ET.Element("MCFile")
        e.set("treename", self.treename)
        e.set("filename", self.filename)
        e.set("scale", str(self.scale))
        e.set("maxevents", str(self.maxevents))
        e.set("pot", str(self.pot))

        for f in self.friend:
            f.to_xml(e)

        for b in self.branch:
            b.to_xml(e)

        return e

@dataclass
class PROXMLMaker:
    mode: str
    model: PROXMLModel
    detector: list[PROXMLDetector] = field(default_factory=list)
    channel: list[PROXMLChannel] = field(default_factory=list)
    plotpot: float = 2e20
    allow_variation_list: list[PROXMLVariation] = field(default_factory=list)
    deny_variation_list: list[PROXMLVariation] = field(default_factory=list)
    mcfile: list[PROXMLMCFile] = field(default_factory=list)

    def to_xml(self):
        elements = []

        e_mode = ET.Element("mode")
        e_mode.set("name", self.mode)
        elements.append(e_mode)

        for d in self.detector:
            elements.append(d.to_xml())

        for c in self.channel:
            elements.append(c.to_xml())

        e_plotpot = ET.Element("plotpot")
        e_plotpot.set("value", str(self.plotpot))
        elements.append(e_plotpot)

        e_variations = ET.Element("variation_list")
        for v in self.allow_variation_list:
            v.to_xml(e_variations, "allow")
        for v in self.deny_variation_list:
            v.to_xml(e_variations, "deny")
        elements.append(e_variations)

        elements.append(self.model.to_xml())

        for m in self.mcfile:
            elements.append(m.to_xml())

        return "\n\n".join([ET.tostring(e, encoding="unicode") for e in elements])
