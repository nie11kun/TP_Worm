import math
import numpy as np

def calculate_grinding_parameters(R1, alfa1, N1, L1, R2, H1, pA1, hN, L21, L22, X0, Z0, B0, C0):
    """
    计算涡轮磨削相关参数
    
    参数：
    R1: float - 涡轮分度圆半径
    alfa1: float - 涡轮加工角度（度）
    N1: int - 涡轮齿数
    L1: float - B轴中心距离砂轮中心的X向距离
    R2: float - 砂轮半径
    H1: float - 砂轮外沿距离分度圆点的高度
    pA1: float - 蜗杆螺旋角度
    hN: int - 蜗杆旋向 (1右旋 -1左旋)
    L21: float - A轴中心距离砂轮中心的Z向距离（正值表示砂轮在A轴回转中心左侧）
    L22: float - B轴中心距离A轴中心（正值表示A轴中心在B轴回转中心左侧）
    X0: float - 磨削中点分度圆点2砂轮对刀点X坐标
    Z0: float - 磨削中点分度圆点2砂轮对刀点Z坐标
    B0: float - 磨削中点分度圆点2砂轮对刀点B轴角度
    C0: float - 磨削中点分度圆点2砂轮对刀点C轴角度
    """
    def find_circle_from_points(p1, p2, p3):
        # 解线性方程组的方法
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3

        # 构造矩阵
        A = np.array([
            [2 * (x2 - x1), 2 * (y2 - y1)],
            [2 * (x3 - x2), 2 * (y3 - y2)]
        ])
        
        B = np.array([
            x2**2 + y2**2 - x1**2 - y1**2,
            x3**2 + y3**2 - x2**2 - y2**2
        ])
        
        # 计算圆心坐标 (a, b)
        center = np.linalg.solve(A, B)
        a, b = center

        # 计算半径 r
        r = np.sqrt((x1 - a)**2 + (y1 - b)**2)

        return (a, b), r
    
    L2 = L21 * math.cos(math.radians(pA1)) + L22  # B轴中心距离砂轮中心的Z向距离（正值表示砂轮在B轴回转中心左侧） (float)
   #  print(f"L2 (B轴中心距离砂轮中心的Z向距离（正值表示砂轮在B轴回转中心左侧）): {L2}")

    # 步骤 1: 计算 L3, R3, 和 beta1
    L3 = L1 + R2 - H1
   #  print(f"L3 (B轴中心到砂轮分度圆点的垂直距离): {L3}")
    R3 = math.sqrt(L3**2 + L2**2)
   #  print(f"R3 (B轴中心和砂轮分度圆的虚拟半径): {R3}")
    beta1 = math.degrees(math.atan(L2 / L3))
   #  print(f"beta1 (B轴中心和砂轮分度圆虚拟半径和垂直距离的夹角): {beta1}")
    
    # 步骤 2: 计算 L4, X1, Z1, X2, Z2
    L4 = 2 * R3 * math.sin(math.radians(alfa1 / 4))
   #  print(f"L4 (磨削起点分度圆点1相对磨削中点分度圆点2的距离): {L4}")
    X1 = L4 * math.cos(math.radians((180 - alfa1 / 2) / 2 - beta1))
    Z1 = L4 * math.sin(math.radians((180 - alfa1 / 2) / 2 - beta1))
   #  print(f"X1 (磨削起点分度圆点1相对磨削中点分度圆点2的X偏移量): {X1}, Z1 (磨削起点分度圆点1相对磨削中点分度圆点2的Z偏移量): {Z1}")
    X2 = X0 - X1
    Z2 = Z0 + Z1
   #  print(f"X2 (磨削起点分度圆点1坐标.X): {X2}, Z2 (磨削起点分度圆点1坐标.Z): {Z2}")
    
    # 步骤 3: 计算 L5, X3, Z3, X4, Z4
    L5 = 2 * R1 * math.sin(math.radians(alfa1 / 4))
   #  print(f"L5 (磨削起点分度圆点2相对磨削起点分度圆点1的距离): {L5}")
    X3 = L5 * math.cos(math.radians((180 - alfa1 / 2) / 2))
    Z3 = L5 * math.sin(math.radians((180 - alfa1 / 2) / 2))
   #  print(f"X3 (磨削起点分度圆点2相对磨削起点分度圆点1的X偏移量): {X3}, Z3 (磨削起点分度圆点2相对磨削起点分度圆点1的Z偏移量): {Z3}")
    X4 = X2 + X3
    Z4 = Z2 - Z3
   #  print(f"X4 (磨削起点分度圆点2坐标.X): {X4}, Z4 (磨削起点分度圆点2坐标.Z): {Z4}")
    
    # 步骤 4: 计算 B 轴和 C 轴角度
    alfa2 = B0 + alfa1 / 2
    beta2 = (C0 + hN * (alfa1 / 2) / (360 / N1) * 360) % 360
   #  print(f"alfa2 (磨削起点B轴角度): {alfa2}, beta2 (磨削起点C轴角度): {beta2}")
    
    # 步骤 5: 计算 L6, X5, Z5, X6, Z6
    L6 = 2 * R3 * math.sin(math.radians(alfa1 / 4))
   #  print(f"L6 (磨削终点分度圆点3相对磨削中点分度圆点2的距离): {L6}")
    if beta1 < alfa1 / 4:
        X5 = L6 * math.sin(math.radians(90 - beta1 - (180 - alfa1 / 2) / 2))
        Z5 = L6 * math.cos(math.radians(90 - beta1 - (180 - alfa1 / 2) / 2))
        X6 = X0 - X5
    else:
        X5 = L6 * math.sin(math.radians(90 - (alfa1 / 2 - beta1) - (180 - alfa1 / 2) / 2))
        Z5 = L6 * math.cos(math.radians(90 - (alfa1 / 2 - beta1) - (180 - alfa1 / 2) / 2))
        X6 = X0 + X5
   #  print(f"X5 (磨削终点分度圆点3相对磨削中点分度圆点2的X偏移量): {X5}, Z5 (磨削终点分度圆点3相对磨削中点分度圆点2的Z偏移量): {Z5}")
    Z6 = Z0 - Z5
   #  print(f"X6 (磨削终点分度圆点3坐标.X): {X6}, Z6 (磨削终点分度圆点3坐标.Z): {Z6}")
    
    # 步骤 6: 计算 L7, X7, Z7, X8, Z8
    L7 = 2 * R1 * math.sin(math.radians(alfa1 / 4))
   #  print(f"L7 (磨削终点分度圆点2相对磨削终点分度圆点3的距离): {L7}")
    X7 = L7 * math.cos(math.radians((180 - alfa1 / 2) / 2))
    Z7 = L7 * math.sin(math.radians((180 - alfa1 / 2) / 2))
   #  print(f"X7 (磨削终点分度圆点2相对磨削终点分度圆点3的X偏移量): {X7}, Z7 (磨削终点分度圆点2相对磨削终点分度圆点3的Z偏移量): {Z7}")
    X8 = X6 + X7
    Z8 = Z6 + Z7
   #  print(f"X8 (磨削终点分度圆点2坐标.X): {X8}, Z8 (磨削终点分度圆点2坐标.Z): {Z8}")
    
    alfa3 = B0 - alfa1 / 2
    beta3 = -1 * hN * alfa1 / (360 / N1) * 360
   #  print(f"alfa3 (磨削终点B轴角度): {alfa3}, beta3 (磨削终点C轴相对起点增量角度): {beta3}")
    
    center, radius = find_circle_from_points((X4, Z4), (X8, Z8), (X0, Z0))

    X9, Z9 = center
    R4 = radius
   #  print(f"X9 (基准圆心X): {X9}, Z9 (基准圆心Z): {Z9}, R4 (基准圆心半径): {R4}")
    
    return X4, Z4, alfa2, beta2, X8, Z8, alfa3, beta3, X9, Z9, R4

# 参数定义
R1 = 110.05              # 涡轮分度圆半径 (float)
alfa1 = 22.3 * 2         # 涡轮加工角度（度） (float)
N1 = 60                  # 涡轮齿数 (int)
L1 = 606.23              # B轴中心距离砂轮中心的X向距离 (float)
R2 = 572.20 / 2          # 砂轮半径 (float)
H1 = 3.72                # 砂轮外沿距离分度圆点的高度 (float)
pA1 = 7.0                # 蜗杆螺旋角度 (float)
hN = 1                   # 蜗杆旋向 (1右旋 -1左旋) (int)
L21 = 75.1328            # A轴中心距离砂轮中心的Z向距离（正值表示砂轮在A轴回转中心左侧） (float)
L22 = -1.145             # B轴中心距离A轴中心（正值表示A轴中心在B轴回转中心左侧） (float)
X0 = 0                   # 磨削中点分度圆点2砂轮对刀点X坐标 (float)
Z0 = 0                   # 磨削中点分度圆点2砂轮对刀点Z坐标 (float)
B0 = 0                   # 磨削中点分度圆点2砂轮对刀点B轴角度 (float)
C0 = 0                   # 磨削中点分度圆点2砂轮对刀点C轴角度 (float)

tG = 3.0                 # 总磨削深度 (float)
fR = 0.02                # 磨削深度步长 (float)

# 计算步数（包含首尾）
num_steps = int(tG / fR)  # 如 5.0/0.02 = 250
# 生成 currentValue 序列，从 5.00 递减到 0.00（保留两位小数）
step_values = [round(tG - i * fR, 2) for i in range(num_steps + 1)]
# 例如： [5.0, 4.98, 4.96, ..., 0.0]

hN_str = "右旋" if hN == 1 else "左旋"

# 打开输出文件
with open("output.txt", "w", encoding="utf-8") as f:
    
    # ----------- 写入 IF 条件语句部分 -----------

    line = f";涡轮分度圆半径: {R1:.4f}mm;\n"
    line += f";涡轮加工角度: {alfa1:.4f}度;\n"
    line += f";涡轮齿数: {N1}齿;\n"
    line += f";蜗杆螺旋角度: {pA1:.4f}度;\n"
    line += f";蜗杆旋向: {hN_str};\n"
    line += f";砂轮直径: {R2*2:.4f}mm;\n"
    line += f";磨削余量: {tG}mm;\n"
    line += ";****************\n"
    line += "DEF REAL currentValue;\n"
    line += "STOPRE;\n"
    line += "currentValue = R210;\n\n"
    f.write(line)

    for i, val in enumerate(step_values):
        # 构造标签，格式：R1_5_00、R1_4_98、...
        label = "CV_" + f"{val:.2f}".replace(".", "_")
        if i == 0:
            # 第一条条件
            line = f"IF (currentValue>={val:.2f}) GOTOF {label};"
        else:
            prev_val = step_values[i - 1]
            line = f"IF (currentValue<{prev_val:.2f}) AND (currentValue>={val:.2f}) GOTOF {label};"
        f.write(line + "\n")
    # 最后一条条件：当 currentValue 小于 0 时
    f.write("IF (currentValue<0) GOTOF CV_0_00;\n\n")
    
    # ----------- 写入每个标签块内容 -----------
    for val in step_values:
        # 标签名：例如 R1_5_00, R1_4_98, ...
        label = "CV_" + f"{val:.2f}".replace(".", "_")

        # 更新当前 R1 值
        current_R1 = R1 - val
        current_X0 = X0 + val

        # 调用计算函数，获得各项磨削参数
        (start_x, start_z, start_b, start_c,
         end_x, end_z, end_b, end_c,
         center_x, center_z, center_r) = calculate_grinding_parameters(
             current_R1, alfa1, N1, L1, R2, H1,
             pA1, hN, L21, L22, current_X0, Z0, B0, C0
         )
        
        # 写入标签及参数块
        f.write(f"{label}:\n")
        f.write(";基准圆心坐标\n")
        f.write(f"R221={center_x:10.4f}        ;X坐标\n")
        f.write(f"R222={center_z:10.4f}        ;Z坐标\n")
        f.write(f"R223={center_r:10.4f}        ;圆弧半径\n")
        f.write(";============\n")
        f.write(";磨削起点坐标\n")
        f.write(f"R231={start_x:10.4f}        ;X坐标\n")
        f.write(f"R232={start_z:10.4f}        ;Z坐标\n")
        f.write(f"R233={start_b:10.4f}        ;B绝对角度\n")
        f.write(f"R234={start_c:10.4f}        ;C绝对角度\n")
        f.write(";============\n")
        f.write(";磨削终点坐标\n")
        f.write(f"R241={end_x:10.4f}        ;X坐标\n")
        f.write(f"R242={end_z:10.4f}        ;Z坐标\n")
        f.write(f"R243={end_b:10.4f}        ;B绝对角度\n")
        f.write(f"R244={end_c:10.4f}        ;C轴终端增量角度\n")
        f.write("RET\n\n")

print("输出已完成，请查看 output.txt 文件。")