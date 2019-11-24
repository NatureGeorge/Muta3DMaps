# @Date:   2019-10-23T22:15:20+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: TableToMD.py
# @Last modified time: 2019-10-23T22:16:44+08:00
import os

tabStr = '''
<div class="output_subarea output_html rendered_html output_result"><div>
<style scoped="">
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Entry</th>
      <th>Gene names</th>
      <th>Status</th>
      <th>Alternative products (isoforms)</th>
      <th>Organism</th>
      <th>Protein names</th>
      <th>canonical_isoform</th>
      <th>unp_map_tage</th>
      <th>yourlist</th>
      <th>UniProt</th>
      <th>GENE</th>
      <th>GENE_status</th>
      <th>Mapping_status</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Q96NU1</td>
      <td>SAMD11</td>
      <td>reviewed</td>
      <td>ALTERNATIVE PRODUCTS:  Event=Alternative promo...</td>
      <td>Homo sapiens (Human)</td>
      <td>Sterile alpha motif domain-containing protein ...</td>
      <td>Q96NU1-3</td>
      <td>Untrusted &amp; No Isoform</td>
      <td>NP_689699</td>
      <td>Q96NU1</td>
      <td>SAMD11</td>
      <td>True</td>
      <td>No</td>
    </tr>
    <tr>
      <th>1</th>
      <td>P43489</td>
      <td>TNFRSF4 TXGP1L</td>
      <td>reviewed</td>
      <td>NaN</td>
      <td>Homo sapiens (Human)</td>
      <td>Tumor necrosis factor receptor superfamily mem...</td>
      <td>NaN</td>
      <td>Trusted &amp; No Isoform</td>
      <td>NP_003318</td>
      <td>P43489</td>
      <td>TNFRSF4</td>
      <td>True</td>
      <td>Yes</td>
    </tr>
    <tr>
      <th>2</th>
      <td>A0A024R084</td>
      <td>SDF4 hCG_19193</td>
      <td>unreviewed</td>
      <td>NaN</td>
      <td>Homo sapiens (Human)</td>
      <td>Stromal cell derived factor 4, isoform CRA_c</td>
      <td>NaN</td>
      <td>Trusted &amp; No Isoform</td>
      <td>NP_057260</td>
      <td>A0A024R084</td>
      <td>SDF4</td>
      <td>True</td>
      <td>No</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Q96L58</td>
      <td>B3GALT6</td>
      <td>reviewed</td>
      <td>NaN</td>
      <td>Homo sapiens (Human)</td>
      <td>Beta-1,3-galactosyltransferase 6 (Beta-1,3-Gal...</td>
      <td>NaN</td>
      <td>Trusted &amp; No Isoform</td>
      <td>NP_542172</td>
      <td>Q96L58</td>
      <td>B3GALT6</td>
      <td>True</td>
      <td>Yes</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Q7RTX0</td>
      <td>TAS1R3 T1R3 TR3</td>
      <td>reviewed</td>
      <td>NaN</td>
      <td>Homo sapiens (Human)</td>
      <td>Taste receptor type 1 member 3 (Sweet taste re...</td>
      <td>NaN</td>
      <td>Trusted &amp; No Isoform</td>
      <td>NP_689414</td>
      <td>Q7RTX0</td>
      <td>TAS1R3</td>
      <td>True</td>
      <td>Yes</td>
    </tr>
  </tbody>
</table>
</div></div>
'''

#
# | Tables        | Are           | Cool  |
# | ------------- |:-------------:| -----:|
# | col 3 is      | right-aligned | $1600 |
# | col 2 is      | centered      |   $12 |
# | zebra stripes | are neat      |    $1 |

divStr = ":-------------:"
maxTdSpace = len(divStr)

# tdStr:td的内容,并且判断中文的个数,因为一个中文和一个字符的len是相同的


def getTdRemainSpaceCount(tdStr):
    charCount = len(tdStr)
    hanziCharCount = 0
    for index in range(0, charCount):
        # 如果是汉字
        if(u'\u4e00' <= tdStr[index] <= u'\u9fff'):
            hanziCharCount += 1
    return maxTdSpace - len(tdStr) - hanziCharCount


def getSpaceStr(spaceCount):
    spaceStr = ""
    for i in range(0, spaceCount):
        spaceStr += " "
    return spaceStr

# 打印一行(一个 <tr>)   targetPreSplitStr 可能为 <td style="text-align: center">, targetBackSplitStr是 </th> 或 </td>


def printTabRow(str, targetPreSplitStr, targetBackSplitStr):
    rawCount = str.count(targetBackSplitStr, 0, len(str))
    splitStr = str
    for index in range(0, rawCount):
        # print("for 循环 splitStr=" + splitStr)
        backIndex = splitStr.find(targetBackSplitStr)
        preStr = splitStr[0:backIndex]
        targetPreIndex = preStr.find(targetPreSplitStr)
        splitedStr = preStr[targetPreIndex:backIndex]
        preIndex = splitedStr.find(">") + 1
        targetStr = splitedStr[preIndex:backIndex]
        if index == 0:
            print("|", end='')

        remainPreSpaceCount = 0
        remainBackSpaceCount = 0
        tempRemainSpaceCount = getTdRemainSpaceCount(targetStr)
        if tempRemainSpaceCount % 2 == 0:
            remainPreSpaceCount = int(tempRemainSpaceCount / 2)
            remainBackSpaceCount = remainPreSpaceCount
        else:
            remainPreSpaceCount = int(tempRemainSpaceCount / 2)
            remainBackSpaceCount = remainPreSpaceCount + 1

        preSpaceStr = getSpaceStr(remainPreSpaceCount)
        backSpaceStr = getSpaceStr(remainBackSpaceCount)
        # 将 strong 标签替换掉
        targetStr = targetStr.replace("<strong>", "**")
        targetStr = targetStr.replace("</strong>", "**")
        print(preSpaceStr + targetStr + backSpaceStr + "|", end='')
        splitStr = splitStr[backIndex + len(targetBackSplitStr):len(str)]
        # print("for 循环，截取 splitStr=" + splitStr)
    print("")  # 如果要是 print("\n")反而会换两行

# divCount 是格数


def printTabDivision(divCount):
    for index in range(0, divCount):
        if index == 0:
            print("|", end='')
        print(divStr + "|", end='')
    print("")


# 解析 thead
tHeadPreIndex = tabStr.find('<thead>') + len('<thead>')
tHeadBackIndex = tabStr.find('</thead>')
tHeadStr = tabStr[tHeadPreIndex:tHeadBackIndex]
# print("tHeadStr=" + tHeadStr)
printTabRow(tHeadStr, "<th", "</th>")
divCount = tHeadStr.count("</th>", 0, len(tHeadStr))
printTabDivision(divCount)

# 解析 tbody
tbodyPreIndex = tabStr.find('<tbody>') + len('<tbody>')
tbodyBackIndex = tabStr.find('</tbody>')
tbodyStr = tabStr[tbodyPreIndex:tbodyBackIndex]
trCount = tbodyStr.count("</tr>", 0, len(tbodyStr))

tempTbodyStr = tbodyStr
for index in range(0, trCount):
    lastTrIndex = tempTbodyStr.find("</tr>") + len("</tr>")
    preTrIndex = tempTbodyStr.find("<tr>")
    targetTrStr = tempTbodyStr[preTrIndex:lastTrIndex]
    printTabRow(targetTrStr, "<td", "</td>")
    tempTbodyStr = tempTbodyStr[lastTrIndex:len(tempTbodyStr)]

# os.system("pause")
