import pandas

# Function to read the mapping file and create a dictionary
def load_taxonomic_mapping(file_path):
        """
        Reads the taxonomic ID to species mapping from a tab-delimited file
        and returns a dictionary with taxonomic IDs as keys and species as values.
        """
        # Read the file using pandas
        df = pandas.read_csv(file_path, sep="\t", header=None, 
                                         names=["code", "kingdom", "taxid", "official", "common", "synonym"])
   
        # Create a dictionary with taxonomic ID as key and species as value
        taxonomic_mapping = dict(zip(df["official"], df["taxid"]))
        return taxonomic_mapping

#function to update a shared dictionary 
def update_shared_dict(shared_dict, key, value):
	shared_copy = shared_dict.get().copy()
	shared_copy[key] = value
	shared_dict.set(shared_copy)

