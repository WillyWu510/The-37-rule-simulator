使用套件如下:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.colors as Colors
from matplotlib import animation
from itertools import combinations
import time

如何重現Demo時的結果:
執行sim_and_particles.py，模擬便開始執行。
程式會分別按10個策略執行100次模擬，每10次模擬計算一次在該策略下的機率。一次模擬約5-20秒
每次模擬結束，程式會print提示信息顯示結果，如下:
--------------------------------------------------------------------------------------------
一個list:顯示在選擇之前遇到那些人
I choose you!!!:做出選擇，並在下一行顯示選擇的人的分數
Yuck, I kiss a real frog.:顯示本次模擬中你並沒有選到王子
Yes! You're my Mr.Right.:顯示本次模擬中你選到王子
No, I throw away the prince. I'm now a frog...Ribbit.:你在你的策略次數之前就碰到王子，所以沒有遇到王子
Too long:模擬執行時間過長(100秒)，之後終止模擬，並將本次結果排除機率之計算。
--------------------------------------------------------------------------------------------
程式大概會跑4-5個小時，並在最後輸出一個list，內容為十個list，每個list的內容為在給定策略下得到的十組機率。

分工表:
吳以理:程式設計、影片拍攝、解讀數據
謝廷鴻:程式設計、影片拍攝及後製、製作ppt
