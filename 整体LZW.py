# -*- coding: gbk -*-
__author__ = 'guolei'
__date__ = '2018/5/21'
# 在上一版基础上将每一层系数LZW改为将各层系数拼在一起整体LZW

def LZW(inStr, narrow=False, bits=14):
    '''使用LZW压缩算法压缩。
        narrow为True时输出字节流位紧缩
        默认最大位宽14位，允许范围12~16位'''
    if isinstance(inStr, str):
        inStr = list(inStr)
        for i in range(len(inStr)):
            inStr[i] = ord(inStr[i])
    sOutStr = [256]  # 先放一个 开始&清除 标记
    mTagMap = {}  # 字典
    iTagCurrent = 258  # 当前的标记总数 0~255 标记256 为begin & clear , 257 为end
    iBitsLen = 9

    iTag = inStr[0]  # 当前标记
    ii = 0
    cTemp = 0
    if bits > 16 or bits < 12:
        return None
    iMaxLen = (1 << bits) - 1  # 最多允许的标记数，由位宽决定

    for ii in range(1, len(inStr)):
        cChar = inStr[ii]
        cTemp = (iTag << 8) + cChar  # （前缀 后缀）
        if cTemp in mTagMap:  # 该（前缀 后缀）已存在
            iTag = mTagMap[cTemp]  # 取出其标记
        else:  # 不存在
            sOutStr.append(iTag)  # 将前缀放入输出流
            mTagMap[cTemp] = iTagCurrent
            iTagCurrent += 1  # 增加一个标记，并放入字典
            iTag = cChar  # 当前标记为后缀

        if iTagCurrent >= iMaxLen:  # 如果到达了最大标记数，清除字典，从头开始
            sOutStr.append(256)
            mTagMap = {}
            iTagCurrent = 258

    if iTag != 0:
        sOutStr.append(iTag)
    sOutStr.append(257)  # 放入结束标记

    if narrow:  # 位紧缩
        return Narrow(sOutStr)
    else:
        return sOutStr

def Narrow(sOutStr):
    sOutN = []
    iTemp = 0
    BitLeft = 0
    nowBit = 9  # 当前位宽
    nowStag = 1 << nowBit  # 当前位宽允许的标记数
    nowTagCount = 258
    for cChar in sOutStr:
        iTemp = iTemp + (cChar << BitLeft)
        nowTagCount += 1
        BitLeft += nowBit

        if cChar == 256:
            nowBit = 9
            nowStag = 1 << nowBit
            nowTagCount = 258

        if nowTagCount >= nowStag:
            nowBit += 1
            nowStag = 1 << nowBit

        while BitLeft >= 8:
            sOutN.append(iTemp & 0xff)
            iTemp = iTemp >> 8
            BitLeft -= 8
    if BitLeft > 0:
        sOutN.append(iTemp)
    return sOutN

def UnNarrow(inStr):
    sOut = []
    iTemp = 0
    BitLeft = 0
    nowBit = 9
    nowStag = 1 << nowBit
    mask = nowStag - 1
    nowTagCount = 258
    for cChar in inStr:
        iTemp = iTemp + (cChar << BitLeft)
        BitLeft += 8
        if BitLeft >= nowBit:
            cTemp = iTemp & mask
            iTemp = iTemp >> nowBit
            BitLeft -= nowBit
            sOut.append(cTemp)
            nowTagCount += 1
            if nowTagCount >= nowStag:
                nowBit += 1
                nowStag = 1 << nowBit
                mask = nowStag - 1
            if cTemp == 256:
                nowBit = 9
                nowStag = 1 << nowBit
                mask = nowStag - 1
                nowTagCount = 258
    if BitLeft > 0:
        sOut.append(iTemp)

    return sOut


def deTag(mTagMap, nowTag, outStr):
    '''''将一个标记转化为元字符序列，并放入到输出流中'''
    if nowTag >= 0:
        sTemp = []
        while nowTag > 255:
            pair = mTagMap[nowTag]
            sTemp.append(pair[1])
            nowTag = pair[0]
        sTemp.append(nowTag)
        sTemp.reverse()
        outStr += sTemp
    return nowTag


def UnLZW(inStr, narrow=False):
    if narrow:
        inStr = UnNarrow(inStr)

    mTagMap = {}
    outStr = []
    nowTagCount = 258
    sTemp = []
    nowTag = -1
    for cChar in inStr:
        if cChar == 256:
            if nowTag >= 0:
                deTag(mTagMap, nowTag, outStr)
            mTagMap = {}
            nowTagCount = 258
            nowTag = -1
        elif cChar == 257:
            if nowTag >= 0:
                deTag(mTagMap, nowTag, outStr)
                nowTag = -1
            return outStr
        elif nowTag == -1:  # 刚开始
            nowTag = cChar
        else:
            pair = [nowTag, 0]
            mTagMap[nowTagCount] = pair
            nowTagCount += 1
            surfix = cChar
            while surfix > 255:
                if surfix not in mTagMap:
                    print('Thera are errors in input')
                    return outStr
                surfix = mTagMap[surfix][0]
            pair[1] = surfix
            deTag(mTagMap, nowTag, outStr)
            nowTag = cChar


def coef_pyramid_plot(coefs, first=0, scale='uniform', ax=None):
    """
    Parameters
    ----------
    coefs : array-like
        Wavelet Coefficients. Expects an iterable in order Cdn, Cdn-1, ...,
        Cd1, Cd0.
    first : int, optional
        The first level to plot.
    scale : str {'uniform', 'level'}, optional
        Scale the coefficients using the same scale or independently by
        level.
    ax : Axes, optional
        Matplotlib Axes instance

    Returns
    -------
    Figure : Matplotlib figure instance
        Either the parent figure of `ax` or a new pyplot.Figure instance if
        `ax` is None.
    """

    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, axisbg='lightgrey')
    else:
        fig = ax.figure

    n_levels = len(coefs)
    n = 2 ** (n_levels - 1)  # assumes periodic

    if scale == 'uniform':
        biggest = [np.max(np.abs(np.hstack(coefs)))] * n_levels
    else:
        # multiply by 2 so the highest bars only take up .5
        biggest = [np.max(np.abs(i)) * 2 for i in coefs]

    for i in range(first, n_levels):
        x = np.linspace(2 ** (n_levels - 2 - i), n - 2 ** (n_levels - 2 - i), 2 ** i)
        ymin = n_levels - i - 1 + first
        yheight = coefs[i] / biggest[i]
        ymax = yheight + ymin
        ax.vlines(x, ymin, ymax, linewidth=1.1)

    ax.set_xlim(0, n)
    ax.set_ylim(first - 1, n_levels)
    ax.yaxis.set_ticks(np.arange(n_levels - 1, first - 1, -1))
    ax.yaxis.set_ticklabels(np.arange(first, n_levels))
    ax.tick_params(top=False, right=False, direction='out', pad=6)
    ax.set_ylabel("Levels", fontsize=14)
    ax.grid(True, alpha=.85, color='white', axis='y', linestyle='-')
    ax.set_title('Wavelet Detail Coefficients', fontsize=16,
                 position=(.5, 1.05))
    fig.subplots_adjust(top=.89)

    return fig

if __name__=='__main__':
    import timeit
    import numpy as np
    import pywt
    import matplotlib.pyplot as plt
    from statsmodels.robust import stand_mad

    # 程序计时开始
    start=timeit.default_timer()

    # 原始数据文件夹20180511_170342.txt
    # 压缩后导出文件夹outFile
    # 位数bit=14
    # 小波基db8
    File='20180511_170342.txt'
    outFile='3层.txt'
    bit=14
    wavtag = 'db8'

    # 原始数据读取
    # 信号显示
    data = []
    with open(r"20180511_170342.txt", "r") as f:
        reader = f.readlines()
        for r in reader:
            r_float = float(r)
            data.append(r_float)
    xx = np.linspace(0, 1, 10240, endpoint=False)
    plt.plot(xx, data)
    plt.title('Original Signal')
    plt.show()

    # 调用wavedec()多级分解函数对数据进行小波变换
    # mode指定了数据补齐的方式
    #‘per’指周期延拓数据
    # 分解层数为9
    wavelet_coefs = pywt.wavedec(data, wavtag, level=3, mode='per')

    # 进行阈值value构造
    # 进行软阈值去噪
    # cA9为细节系数
    # cD9, cD8, cD7, cD6, cD5, cD4, cD3, cD2, cD1 为模糊系数
    # denoised_coefs格式为[array([...]),array([...]),...,array([...])]
    sigma = stand_mad(wavelet_coefs[-1])
    uthresh = sigma * np.sqrt(2 * np.log(len(data)))
    denoised_coefs = wavelet_coefs[:]
    denoised_coefs[1:] = (pywt.threshold(coefs, value=uthresh, mode='soft') for coefs in denoised_coefs[1:])
    #cA9, cD9, cD8, cD7, cD6, cD5, cD4, cD3, cD2, cD1 = denoised_coefs

    # 将去噪得到的系数进行处理，以方便存入txt
    # 首先将denoised_coefs中每个array转换tolist为列表，在append入noisy_coefs_list
    # 再将noisy_coefs_list中每个列表的float数据转换为str并以“ ”分隔
    # 将i_str进行LZW压缩
    # 返回的aIn为list，将其中每个0~256整数做参数，返回一个对应字符并拼接
    wavelet_coefs_list=[]
    for i in denoised_coefs:# i为array
        wavelet_coefs_list.append(i.tolist())
    storeData=''
    i_strAll=''
    for i in wavelet_coefs_list:
        i_str=" ".join(map(str,i))
        i_str=i_str+"!"
        i_strAll+=i_str


    dataLZW = LZW(i_strAll, False)
    for j in dataLZW:
        storeData+=chr(j)

    # 打开outFile
    # 将storeData先编码后写入
    f = open(outFile, 'wb')
    try:
        f.write(storeData.encode())
    finally:
        f.close()

    # 程序计时结束
    end = timeit.default_timer()

    # 打开outFile
    # 读取已经压缩的数据进行下一步解压缩
    with open(outFile,'rb') as f:
        readerOutFile=f.read().decode()

    # 转换回数字并加入列表a
    # 以257（结束标志）为分隔将一组系数加入列表signalData
    # 将signalData传入LZW解压函数进行解压
    # 返回的bData为str，将所有系数放入列表signalDataFin
    dataEncode = []
    for i in readerOutFile:
        dataEncode.append(ord(i))
    unLZWDataNumb=UnLZW(dataEncode, False)
    dataFin = ''
    for i in unLZWDataNumb:
        dataFin += chr(i)
    dataTemp=''
    signalFin = []
    for i in dataFin:
        if i!="!":
            dataTemp+=i
        else:
            listTempStr = dataTemp.split(" ")
            listTempFloat = [float(i) for i in listTempStr]
            listTempArray = np.array(listTempFloat)
            signalFin.append(listTempArray)
            dataTemp=''

    # 进行信号重建
    signal = pywt.waverec(signalFin, wavtag, mode='per')

    # 画原始信号图像和重建图像对比
    fig = plt.figure(7)
    ax71 = fig.add_subplot(211)
    ax71.plot(xx, data)
    ax71.grid(True)
    ax71.set_title('Original Signal')
    ax71.tick_params(labelbottom=False)

    ax72 = fig.add_subplot(212)
    ax72.plot(xx, signal)
    ax72.grid(True)
    ax72.set_title('Denoised Signal')
    ax72.legend()

    plt.show()

    # 输出程序运行时间
    print(end-start)

    # 程序结束
    print('end')
