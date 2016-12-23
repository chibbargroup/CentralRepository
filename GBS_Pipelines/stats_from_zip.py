import zipfile,os.path, math, sys
#This script will generate useful information in a spreadsheet when run on a fastqc zip file


def unzip(source_filename, dest_dir):
    with zipfile.ZipFile(source_filename) as zf:
        for member in zf.infolist():
            # Path traversal defense copied from
            # http://hg.python.org/cpython/file/tip/Lib/http/server.py#l789
            words = member.filename.split('/')
            path = dest_dir
            for word in words[:-1]:
                drive, word = os.path.splitdrive(word)
                head, word = os.path.split(word)
                if word in (os.curdir, os.pardir, ''): continue
                path = os.path.join(path, word)
            zf.extract(member, path)

def delete_dir(top):
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(top)

def find_data_file(zip_root):
    for root, dirs, files in os.walk(zip_root, topdown=False):
        for name in files:
            if name == "fastqc_data.txt":
                return os.path.join(root, name)

def find_average_length(file_path):
    found = False
    totalVal = 0
    totalCount = 0;
    for line in open(file_path, 'r'):
        if ">>END_MODULE" in line:
            found = False
        if found and not "#" in line:
            values = line.rstrip().split("\t")
            minLength = values[0].split("-")[0]
            maxLength = values[0].split("-")[1]
            avgLength = (float(maxLength) + float(minLength)) / 2.0
            numElements = float(values[1])
            totalCount = totalCount + numElements
            totalVal = totalVal + (avgLength * numElements)
        if ">>Sequence Length Distribution" in line:
            found = True
    totalAvg = totalVal / totalCount
    return totalAvg

def find_std_dev_length(file_path, average):
    found = False
    sumSquared = 0
    totalElements = 0
    for line in open(file_path, 'r'):
        if ">>END_MODULE" in line:
            found = False
        if found and not "#" in line:
            values = line.rstrip().split("\t")
            minLength = values[0].split("-")[0]
            maxLength = values[0].split("-")[1]
            avgLength = (float(maxLength) + float(minLength)) / 2.0
            numElements = float(values[1])
            diff = average - avgLength
            diffSquared = diff * diff
            sumSquared = sumSquared + diffSquared * numElements
            totalElements = totalElements + numElements
        if ">>Sequence Length Distribution" in line:
            found = True
    variance = sumSquared / (totalElements - 1)
    stdDev = math.sqrt(variance)
    return stdDev

def find_average_quality(file_path):
    found = False
    totalVal = 0
    totalCount = 0;
    for line in open(file_path, 'r'):
        if ">>END_MODULE" in line:
            found = False
        if found and not "#" in line:
            values = line.rstrip().split("\t")
            quality = float(values[0])
            num = float(values[1])
            totalVal = totalVal + quality * num;
            totalCount = totalCount + num
        if ">>Per sequence quality scores" in line:
            found = True
    totalAvg = totalVal / totalCount
    return totalAvg

def find_std_dev_quality(file_path, average):
    found = False
    sumSquared = 0
    totalElements = 0
    for line in open(file_path, 'r'):
        if ">>END_MODULE" in line:
            found = False
        if found and not "#" in line:
            values = line.rstrip().split("\t")
            quality = float(values[0])
            num = float(values[1])
            diff = average - quality
            diffSquared = diff * diff
            sumSquared = sumSquared + diffSquared * num
            totalElements = totalElements + num
        if ">>Per sequence quality scores" in line:
            found = True
    variance = sumSquared / (totalElements - 1)
    stdDev = math.sqrt(variance)
    return stdDev

def find_total_reads(file_path):
    found = False
    for line in open(file_path, 'r'):
        if ">>END_MODULE" in line:
            found = False
        if found and not "#" in line:
            values = line.rstrip().split("\t")
            if values[0] == "Total Sequences":
                return float(values[1])
        if ">>Basic Statistics" in line:
            found = True
    return  float(-1)

def find_length_range(file_path):
    found = False
    for line in open(file_path, 'r'):
        if ">>END_MODULE" in line:
            found = False
        if found and not "#" in line:
            values = line.rstrip().split("\t")
            if values[0] == "Sequence length":
                return values[1]
        if ">>Basic Statistics" in line:
            found = True
    return  "-"

def find_stats(file_path):
    avgLength = find_average_length(data)
    stdLength = find_std_dev_length(data, avgLength)
    lengthRange = find_length_range(data)
    avgQuality = find_average_quality(data)
    stdQuality = find_std_dev_quality(data, avgQuality)
    totalReads = find_total_reads(data)
    return {"Average Length":avgLength, "Standard Deviation of Length":stdLength, "Range of Length":lengthRange, "Average Quality":avgQuality, "Standard Deviation of Quality":stdQuality, "Total Read Count":totalReads}


if len(sys.argv) > 1:
    rootDir = sys.argv[1]
else:
    rootDir = "./"

statsDict = {}
for root, dirs, files in os.walk(rootDir):
    for file in files:
        filename, file_extension = os.path.splitext(file)
        if file_extension == ".zip" and "_fastqc" in filename:
            filepath = os.path.join(root, file)
            temporaryPath = os.path.join(root, "tmp")
            unzip(filepath, temporaryPath)
            data = find_data_file(temporaryPath)
            stats = find_stats(data)
            statsDict[filename] = stats
            delete_dir(temporaryPath)
#Print Results in tsv format
order = ["Total Read Count","Average Length", "Standard Deviation of Length", "Range of Length", "Average Quality", "Standard Deviation of Quality"]
title_row = "FileID"
for colID in order:
    title_row = title_row + "\t"+colID
print title_row
for fileID in statsDict:
    row_string = fileID
    stats = statsDict[fileID]
    for colID in order:
        row_string = row_string + "\t" + str(stats[colID])
    print row_string