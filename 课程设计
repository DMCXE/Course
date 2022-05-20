from iapws import IAPWS97
import numpy as np
class CoolentWater(object):
    Tin = 25                #入口，摄氏度,20-25
    Tout = 45               #出口，摄氏度,40-50
    Pressure = 0.2          #压力，MPa,0.1-0.2
    Tref = 0.5*(Tin+Tout)   #定性温度（摄氏度）
    water = IAPWS97(T=Tref+273.15,P=Pressure)
    rho = water.Liquid.rho  #密度，kg/m3
    cp = water.Liquid.cp    #定压比热容cp2,kJ/kg·K
    k = water.Liquid.k      #传热系数,W/m·K
    mu = water.Liquid.mu    #动力粘度,Pa·s,kg/(m`s)
    Pr = 1000*cp*mu/k       #普朗特数

class HotWater(object):
    Tin = 80                #入口，摄氏度,75-85
    Tout = 40               #出口，摄氏度,40-50
    Pressure = 0.15         #压力，MPa, 0.1-0.3
    G = 18000               #热水流量，kg/h
    Tref = (Tin+Tout)*0.5   #定性温度（摄氏度）
    water = IAPWS97(T=Tref+273.15,P=Pressure)
    rho = water.Liquid.rho  # 密度，kg/m3
    cp = water.Liquid.cp  # 定压比热容cp1,kJ/kg·K
    k = water.Liquid.k  # 传热系数,W/m·K
    mu = water.Liquid.mu  # 动力粘度,Pa·s,kg/(m`s)
    Pr = 1000 * cp * mu / k  # 普朗特数
'''换热器效率'''
TransferEfficiency =((CoolentWater.Tout-CoolentWater.Tin)
                     if (CoolentWater.Tout-CoolentWater.Tin)>(HotWater.Tin-HotWater.Tout)
                     else HotWater.Tin-HotWater.Tout)/(HotWater.Tin-CoolentWater.Tin)
'''设计传热量Q0/W'''
TransferQ = HotWater.G * HotWater.cp * (HotWater.Tin-HotWater.Tout)*TransferEfficiency*1000/3600
'''冷却水流量G2/Kg/h'''
CoolentWaterQ = (3600*TransferQ)/(1000*CoolentWater.cp*(CoolentWater.Tout-CoolentWater.Tin))
'''逆流平均温差tN/摄氏度'''
Tb = (CoolentWater.Tout-CoolentWater.Tin) \
    if (CoolentWater.Tout-CoolentWater.Tin)>(HotWater.Tin-HotWater.Tout) \
    else (HotWater.Tin-HotWater.Tout)
Ts = (CoolentWater.Tout-CoolentWater.Tin) \
    if (CoolentWater.Tout-CoolentWater.Tin)<(HotWater.Tin-HotWater.Tout) \
    else (HotWater.Tin-HotWater.Tout)
ContraFlowDT = (Tb-Ts)/np.log(Tb/Ts)
'''24.温差校正系数P、R、phi'''
P = (CoolentWater.Tout-CoolentWater.Tin)/(HotWater.Tin-CoolentWater.Tin)
R = (HotWater.Tin-HotWater.Tout)/(CoolentWater.Tout-CoolentWater.Tin)
phi = 0.9
'''25.有效平均温差'''
EffectiveT = phi*ContraFlowDT
'''26-34采用四舍五入取整'''
'''26.试选传热系数K0'''
K0 = 500
'''27.初选传热面积F0'''
F0 = TransferQ/(K0*EffectiveT)
'''28.管子外径d0,内径di'''
d0 = 25*10**(-3)
di = d0-2*2.5*10**(-3)
'''29.管子长度l'''
pipelength = 2.5
'''30.总管数Nt'''
TotPipe = round(F0/(np.pi*d0*pipelength)) #向上取整
'''31.管程流通横截面a2'''
PipeFlowSection = (0.5*TotPipe)*0.25*np.pi*di*di

'''32.管程流速w2'''
w2 = CoolentWaterQ/(CoolentWater.rho*PipeFlowSection*3600)
'''33.管程雷诺数Re2'''
Re2 = CoolentWater.rho*w2*di/CoolentWater.mu
'''34.管程换热系数h2'''
Nu2 = CoolentWater.rho*w2/CoolentWater.mu
h2 = (0.023*Nu2**0.8*CoolentWater.Pr**0.3)*CoolentWater.k/di
'''35.管子排列方式'''
#正三角形
'''36.管间距S/m'''
S = 32*10**(-3)
'''37.管束中心处一排管束Nc'''
Nc = round(1.1*np.sqrt(TotPipe))
'''38.管束外沿与壳体间距'''
e = 2*d0
'''39.壳体内径'''
Ds = S*(Nc-1)+4*d0
'''40.长径比'''
#参考值4-25，最佳值6-10
LRrate = pipelength/Ds
'''41.弓形折流板弓高h'''
ArcH = 0.2*Ds
'''42.折流板间距B'''
ArcB = Ds/3
'''43.折流板数nb'''
ArcN = round((pipelength/ArcB)-1)
'''44.壳程流通截面a1'''
ShellSection = ArcB*Ds*(1-d0/S)
'''45.壳程流速w1'''
ShellV = HotWater.G/(3600*HotWater.rho*ShellSection)
'''46.壳程质量流速W1'''
ShellG = HotWater.rho*ShellV
'''47.壳程当量直径de'''
ShellDe = (3.464*S**2-np.pi*d0**2)/(np.pi*d0)
'''48.壳程雷诺数Rel'''
ShellRel = ShellV*ShellDe/HotWater.mu
'''49.管间距比值'''
PipDisRate = 2/np.sqrt(3)

'''计算排数'''
PipP = round((4*Nc-1)/5)
'''51.管排修正系数'''
Piptable = [0.64,0.76,0.84,0.89,0.93,0.945,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.993,0.996]
if PipP<16:
    Pipk = Piptable[PipP-1]
else: Pipk = 1
'''50.努塞尔数Nu1/不考虑物性变化'''
if PipP>=16:
    if ShellRel>1 and ShellRel<500:
        ShellNuf = 1.04*(ShellRel)**0.4*HotWater.Pr**0.36
    if ShellRel>500 and ShellRel<1000:
        ShellNuf = 0.71 * (ShellRel) ** 0.5 * HotWater.Pr ** 0.36
    if ShellRel>1000 and ShellRel<2*10**5:
        ShellNuf = 0.35*PipDisRate**0.2 * (ShellRel) ** 0.6 * HotWater.Pr ** 0.36
    if ShellRel>2*10**5 and ShellRel< 2*10**6:
        ShellNuf = 0.031*PipDisRate**0.2 * (ShellRel) ** 0.8 * HotWater.Pr ** 0.36
else:
    if ShellRel>1 and ShellRel<500:
        ShellNuf = Pipk*1.04*(ShellRel)**0.4*HotWater.Pr**0.36
    if ShellRel>500 and ShellRel<1000:
        ShellNuf = Pipk* 0.71 * (ShellRel) ** 0.5 * HotWater.Pr ** 0.36
    if ShellRel>1000 and ShellRel<2*10**5:
        ShellNuf = Pipk* 0.35*PipDisRate**0.2 * (ShellRel) ** 0.6 * HotWater.Pr ** 0.36
    if ShellRel>2*10**5 :#and ShellRel< 2*10**6:
        ShellNuf = Pipk* 0.031*PipDisRate**0.2 * (ShellRel) ** 0.8 * HotWater.Pr ** 0.36

'''52.壳程换热系数h1'''
Shellh = Pipk * ShellNuf * HotWater.k /d0
'''53.冷却水侧污垢热阻r2/以蒸馏水为例'''
CoolentWaterR = 8.6*10**(-5)
'''54.热水侧污垢热阻r1/以蒸馏水为例'''
HotWaterR = 8.6*10**(-5)
'''56.总污垢热阻r1/以蒸馏水为例'''
Rsigma = Shellh**(-1)+HotWaterR+CoolentWaterR*(d0/di)+Pipk**(-1)*(d0/di)
'''57/传热系数'''
Kj = 1/Rsigma
'''58.传热系数比值大小,1.1-1.2'''
KRate = Kj/K0
'''59.管外壁热流密度ql'''
ShellQ = TransferQ/(TotPipe*np.pi*d0*pipelength)
'''60.管外壁温度tw1'''
ShellT = HotWater.Tref-ShellQ*(1/Shellh+HotWaterR)
print(KRate)
print(PipeFlowSection)
