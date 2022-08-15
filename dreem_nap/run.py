import sys, os

path = os.path.dirname('/'.join(os.path.abspath(__file__).split('/')[:-1]))
sys.path.append(path)

def main(**args):
    """
    DREEM processes DMS next generation sequencing data to produce mutational
    profiles that relate to DMS modification rates written by Silvi Rouskin and the
    Rouskin lab (https://www.rouskinlab.com/)
    """
    run(args)

def run(args):
    print('Done!')

if __name__ == "__main__":
    sys.argv = sys.argv[0]
    main()