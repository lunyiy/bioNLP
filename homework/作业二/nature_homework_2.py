#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: zhiyan

import matplotlib.pylab as plt
from wordcloud import WordCloud
import jieba
# import numpy as np
# import PIL .Image as image

# 打开词云文件，获取内容
f = open('e:/VScode/python/total-GENiA.txt', "r", encoding="utf-8")
text = f.read()
f.close()
splitWord = jieba.cut(text)

# 获取特定形状的图片
# mask = np.array(image.open("e:/VScode/python/template1.jpg"))

wordcloud = WordCloud(
    background_color="white",
    width=800,
    height=550
    # mask=mask
    ).generate(' '.join(splitWord))

plt.imshow(wordcloud, interpolation='bilinear')
plt.axis("off")
plt.show()