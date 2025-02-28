import matplotlib.pyplot as plt
import csv
from dna_features_viewer import BiopythonTranslator
# Script modified from custom_biopython_translator.py
# and with_gc_plot.py


class CustomTranslatorWithNewFeatures(BiopythonTranslator):
    """Custom translator implementing the following theme:

    - Color terminators in green, CDS in blue, all other features in gold.
    - Do not display features that are restriction sites unless they are BamHI
    - Do not display labels for restriction sites.
    - For CDS labels just write "CDS here" instead of the name of the gene.
    - Also draw new features called "origin", "model_pause" and "fork_stall".
    """

    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            return "blue"
        elif feature.type == "terminator":
            return "green"
        else:
            return "gold"

    def compute_feature_label(self, feature):
        if feature.type == 'restriction_site' or feature.type == 'misc_feature':
            return None
        elif feature.type == "CDS":
            return "CDS here"
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

    def compute_filtered_features(self, features):
        """Do not display promoters. Just because."""
        return [
            feature for feature in features
            if (feature.type != "restriction_site")
               or ("BamHI" in str(feature.qualifiers.get("label", '')))
        ]

    def translate_feature(self, feature):
        """ Add new features called "origin", "model_pause" and "fork_stall". """

        # Translate the feature as usual
        t_feature = super().translate_feature(feature)

        # now add the new features
        if feature.type == "origin":
            t_feature.data['type'] = "origin"
            t_feature.color = "yellow"
            t_feature.thickness = 0.4
            t_feature.start = t_feature.start - 20
            t_feature.end = t_feature.end + 20
        elif feature.type == "model_pause":
            t_feature.data['type'] = "model_pause"
            t_feature.color = "brown"
            t_feature.thickness = 0.4
            t_feature.label = None
            t_feature.start = t_feature.start - 20
            t_feature.end = t_feature.end + 20
        elif feature.type == "fork_stall":
            t_feature.data['type'] = "fork_stall"
            t_feature.color = "green"
            t_feature.thickness = 0.4
            t_feature.label = None
            t_feature.start = t_feature.start - 40
            t_feature.end = t_feature.end + 40

        # show how to draw a strandless feature without a label
        if t_feature.label == "kanR terminator":
            t_feature.strand = 0
            t_feature.label = None

        return t_feature


fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, figsize=(10, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1, 1]}
)

graphic_record = CustomTranslatorWithNewFeatures().translate_record("example_sequence_new_features.gb")
graphic_record.crop((2000, 4000)).plot(figure_width=10, ax=ax1, with_ruler=False)

x2_bedgraph = []
y2_bedgraph = []
with open('bedgraph_a', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='#')
    for row in reader:
        x2_bedgraph.append(int(row[1])/2 + int(row[2])/2)
        y2_bedgraph.append(float(row[3]))

ax2.fill_between(x2_bedgraph, y2_bedgraph, alpha=0.3)
ax2.set_ylim(bottom=0)
ax2.set_ylabel("signal_1")
ax2.set_xlabel("Position (bp)")

x3_bedgraph = []
y3_bedgraph = []
with open('bedgraph_b', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='#')
    for row in reader:
        x3_bedgraph.append(int(row[1])/2 + int(row[2])/2)
        y3_bedgraph.append(float(row[3]))

ax3.fill_between(x3_bedgraph, y3_bedgraph, alpha=0.3)
ax3.set_ylim(bottom=50)
ax3.set_ylabel("signal_2")
ax3.set_xlabel("Position (bp)")

fig.tight_layout()  # Resize the figure to the right height
fig.savefig("new_features_with_bedgraph.png")
