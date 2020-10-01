import os

def main():
    path = "./data/"
    coord_files = os.listdir(path)
    
    print()
    
    for coord_file in coord_files:
        d = count_fragments_per_entry(path + coord_file)
        print(coord_file)
        print("Entries: ", len(d))
        print("All fragments: ", sum(d.values()))

        entries_more_than_one = 0
        for value in d.values():
            if value > 1:
                entries_more_than_one += 1

        print("entries with more than one fragment: ", entries_more_than_one, end= "\n\n")



def count_fragments_per_entry(filename):
    with open(filename) as inputfile:
        lines = inputfile.readlines()
    
    count_dict = {}

    for line in lines:
        if "FRAG" in line:
            information = line.split('**')
            entry = information[0].strip()

            if entry in count_dict.keys():
                count_dict[entry] += 1
            else:
                count_dict[entry] = 1

    return count_dict

if __name__ == "__main__":
    main()
