# Splits the list of strains into SEGMENT_COUNT number of files

import sys

SEGMENT_COUNT = 4

def splitSampleFile(fileName, segmentCount):
    with open(f'{fileName}.txt', 'r') as f: 
        lineCount = 0
        for line in f:
            lineCount += 1

    linesPerSegment = lineCount // segmentCount

    segments = [[] for _ in range(segmentCount)]
    with open(f'{fileName}.txt', 'r') as f:
        segment = 0
        for lineNum, line in enumerate(f):
            segments[segment].append(line.strip())
            if (lineNum + 1) % linesPerSegment == 0:
                if lineNum + linesPerSegment <= lineCount:
                    segment += 1

    for segmentNum, segment in enumerate(segments):
        with open(f'{fileName} [{segmentNum}].txt', 'w') as f:
            for sample in segment:
                f.write(f'{sample}\n')

if __name__ == '__main__':
    #- splitSampleFile(sys.argv[1], sys.argv[2])
    splitSampleFile('./Resources/samples', 4)
