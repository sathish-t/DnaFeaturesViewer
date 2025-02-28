from dna_features_viewer import BiopythonTranslator
# Script modified from custom_biopython_translator.py


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
            t_feature.start = t_feature.start - 80
            t_feature.end = t_feature.end + 80
        elif feature.type == "model_pause":
            t_feature.data['type'] = "model_pause"
            t_feature.color = "brown"
            t_feature.thickness = 0.4
            t_feature.label = None
            t_feature.start = t_feature.start - 80
            t_feature.end = t_feature.end + 80
        elif feature.type == "fork_stall":
            t_feature.data['type'] = "fork_stall"
            t_feature.color = "green"
            t_feature.thickness = 0.4
            t_feature.label = None
            t_feature.start = t_feature.start - 40
            t_feature.end = t_feature.end + 40
        return t_feature


graphic_record = CustomTranslatorWithNewFeatures().translate_record("example_sequence_new_features.gb")
ax, _ = graphic_record.plot(figure_width=10)
ax.figure.tight_layout()
ax.figure.savefig("new_features.png")
