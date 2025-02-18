import os, shutil, platform, subprocess
import math
import numpy as np
from datetime import datetime

def find_circle_from_points(p1, p2, p3):
    """
    已知平面上三个点坐标，求过此三点的圆心和半径。
    
    参数：
        p1, p2, p3: (float, float) - 平面上三点 (x, y)
    返回：
        center: (float, float) - 圆心 (x, y)
        radius: float          - 圆半径
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    A = np.array([
        [2 * (x2 - x1), 2 * (y2 - y1)],
        [2 * (x3 - x2), 2 * (y3 - y2)]
    ])
    
    B = np.array([
        x2**2 + y2**2 - x1**2 - y1**2,
        x3**2 + y3**2 - x2**2 - y2**2
    ])
    
    center = np.linalg.solve(A, B)  # 求圆心 (a, b)
    a, b = center
    radius = math.sqrt((x1 - a) ** 2 + (y1 - b) ** 2)  # 求半径

    return (a, b), radius

def calculate_grinding_parameters(
    R1, alfa1, N1, L1, R2, H1, pA1, hN,
    L21, L22, X0, Z0, B0, C0
):
    """
    计算涡轮磨削相关参数，并返回磨削起点、磨削终点、基准圆心等坐标和角度。
    
    参数：
        R1:  涡轮分度圆半径
        alfa1: 涡轮加工角度（度）
        N1:  涡轮齿数
        L1:  B轴中心距离砂轮中心的 X 向距离
        R2:  砂轮半径
        H1:  砂轮外沿距离分度圆点的高度
        pA1: 蜗杆螺旋角度
        hN:  蜗杆旋向 (1 表示右旋，-1 表示左旋)
        L21: A 轴中心距离砂轮中心的 Z 向距离
        L22: B 轴中心距离 A 轴中心
        X0:  磨削中点分度圆点 2 砂轮对刀点 X 坐标
        Z0:  磨削中点分度圆点 2 砂轮对刀点 Z 坐标
        B0:  磨削中点分度圆点 2 砂轮对刀点 B 轴角度
        C0:  磨削中点分度圆点 2 砂轮对刀点 C 轴角度

    返回：
        (start_x, start_z, start_b, start_c,
         end_x, end_z,   end_b,   end_c,
         center_x, center_z, center_r)
    """
    # 计算 B 轴中心与砂轮中心在 Z 向上的合成距离 L2（正值表示砂轮在B轴回转中心左侧）
    L2 = L21 * math.cos(math.radians(pA1)) + L22

    # 第 1 步：计算 L3, R3, beta1
    L3 = L1 + R2 - H1
   #  print(f"L3 (B轴中心到砂轮分度圆点的垂直距离): {L3}")
    R3 = math.sqrt(L3**2 + L2**2)
   #  print(f"R3 (B轴中心和砂轮分度圆的虚拟半径): {R3}")
    beta1 = math.degrees(math.atan(L2 / L3))
   #  print(f"beta1 (B轴中心和砂轮分度圆虚拟半径和垂直距离的夹角): {beta1}")

    # 第 2 步：磨削中点到磨削起点（在分度圆内）的位移
    L4 = 2 * R3 * math.sin(math.radians(alfa1 / 4))
   #  print(f"L4 (磨削起点分度圆点1相对磨削中点分度圆点2的距离): {L4}")
    X1 = L4 * math.cos(math.radians((180 - alfa1 / 2) / 2 - beta1))
    Z1 = L4 * math.sin(math.radians((180 - alfa1 / 2) / 2 - beta1))
   #  print(f"X1 (磨削起点分度圆点1相对磨削中点分度圆点2的X偏移量): {X1}, Z1 (磨削起点分度圆点1相对磨削中点分度圆点2的Z偏移量): {Z1}")
    X2 = X0 - X1
    Z2 = Z0 + Z1
   #  print(f"X2 (磨削起点分度圆点1坐标.X): {X2}, Z2 (磨削起点分度圆点1坐标.Z): {Z2}")

    # 第 3 步：分度圆半径 R1 决定的磨削起点再向内侧的位移
    L5 = 2 * R1 * math.sin(math.radians(alfa1 / 4))
   #  print(f"L5 (磨削起点分度圆点2相对磨削起点分度圆点1的距离): {L5}")
    X3 = L5 * math.cos(math.radians((180 - alfa1 / 2) / 2))
    Z3 = L5 * math.sin(math.radians((180 - alfa1 / 2) / 2))
   #  print(f"X3 (磨削起点分度圆点2相对磨削起点分度圆点1的X偏移量): {X3}, Z3 (磨削起点分度圆点2相对磨削起点分度圆点1的Z偏移量): {Z3}")
    X4 = X2 + X3
    Z4 = Z2 - Z3
   #  print(f"X4 (磨削起点分度圆点2坐标.X): {X4}, Z4 (磨削起点分度圆点2坐标.Z): {Z4}")

    # 第 4 步：起点 B 轴和 C 轴绝对角度
    alfa2 = B0 + alfa1 / 2
    # 这里的 360/(N1) 对应一齿所占的角度；乘以 (alfa1/2) 再乘以旋向
    beta2 = (C0 + hN * (alfa1 / 2) / (360 / N1) * 360) % 360
   #  print(f"alfa2 (磨削起点B轴角度): {alfa2}, beta2 (磨削起点C轴角度): {beta2}")

    # 第 5 步：从磨削中点到磨削终点（分度圆外侧）的位移
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

    # 第 6 步：分度圆半径 R1 决定的磨削终点再向外侧的位移
    L7 = 2 * R1 * math.sin(math.radians(alfa1 / 4))
   #  print(f"L7 (磨削终点分度圆点2相对磨削终点分度圆点3的距离): {L7}")
    X7 = L7 * math.cos(math.radians((180 - alfa1 / 2) / 2))
    Z7 = L7 * math.sin(math.radians((180 - alfa1 / 2) / 2))
   #  print(f"X7 (磨削终点分度圆点2相对磨削终点分度圆点3的X偏移量): {X7}, Z7 (磨削终点分度圆点2相对磨削终点分度圆点3的Z偏移量): {Z7}")
    X8 = X6 + X7
    Z8 = Z6 + Z7
   #  print(f"X8 (磨削终点分度圆点2坐标.X): {X8}, Z8 (磨削终点分度圆点2坐标.Z): {Z8}")

    # 终点 B 轴角度与 C 轴角度增量
    alfa3 = B0 - alfa1 / 2
    beta3 = -1 * hN * alfa1 / (360 / N1) * 360
   #  print(f"alfa3 (磨削终点B轴角度): {alfa3}, beta3 (磨削终点C轴相对起点增量角度): {beta3}")

    # 计算基准圆心 (X9, Z9) 及其半径 R4
    (X9, Z9), R4 = find_circle_from_points((X4, Z4), (X8, Z8), (X0, Z0))
   #  print(f"X9 (基准圆心X): {X9}, Z9 (基准圆心Z): {Z9}, R4 (基准圆心半径): {R4}")

    # 返回需要的磨削位置和圆心信息
    return (
        X4, Z4, alfa2, beta2,    # 磨削起点
        X8, Z8, alfa3, beta3,    # 磨削终点
        X9, Z9, R4               # 基准圆心坐标和半径
    )

def generate_step_values(total_depth, step_size):
    """
    生成从 total_depth （向上取整到 step_size 整倍数）递减到 0 的磨削步进序列。
    每一步相差 step_size，返回值保留两位小数。
    
    参数：
        total_depth: float - 磨削总深度（若非 step_size 的整倍数，将被适当放大）
        step_size:   float - 磨削步长
    返回：
        values: list[float]
    """
    # 先根据 step_size 将 total_depth 向上取整到整倍数
    # 例如 total_depth=4.95, step_size=0.02，则放大到4.96
    total_depth_rounded = math.ceil(total_depth / step_size) * step_size
    
    # 计算步数（含首尾）
    num_steps = int(total_depth_rounded / step_size)
    
    # 生成从 total_depth_rounded 递减到 0 的序列
    values = [round(total_depth_rounded - i * step_size, 2) for i in range(num_steps + 1)]
    return values

def write_output_file(filename, step_values, param_dict):
    """
    将生成的磨削指令写入文件。
    
    参数：
        filename: str            - 输出文件名
        step_values: List[float] - 磨削深度序列
        param_dict: dict         - 存放各种参数的字典
    """
    # 从字典中读取参数
    (R1, alfa1, N1, L1, R2, H1,
     pA1, hN, L21, L22,
     X0, Z0, B0, C0,
     tG, fR, Ypos) = (
       param_dict["R1"],   param_dict["alfa1"], param_dict["N1"],
       param_dict["L1"],   param_dict["R2"],    param_dict["H1"],
       param_dict["pA1"],  param_dict["hN"],
       param_dict["L21"],  param_dict["L22"],
       param_dict["X0"],   param_dict["Z0"],
       param_dict["B0"],   param_dict["C0"],
       param_dict["tG"],   param_dict["fR"],
       param_dict["Ypos"]
    )

    hN_str = "右旋" if hN == 1 else "左旋"

    # 获取当前日期时间
    current_datetime = datetime.now()

    # 格式化为字符串
    current_datetime_str = current_datetime.strftime("%Y-%m-%d %H:%M:%S")

    with open(filename, "w", encoding="utf-8") as f:
        # 前置信息和 IF 条件
        header = (
            f";涡轮分度圆半径: {R1:.4f}mm;\n"
            f";涡轮加工角度: {alfa1:.4f}度;\n"
            f";涡轮齿数: {N1}齿;\n"
            f";蜗杆螺旋角度: {pA1:.4f}度;\n"
            f";蜗杆旋向: {hN_str};\n"
            f";砂轮直径: {R2*2:.4f}mm;\n"
            f";分度圆齿顶高: {H1:.4f}mm;\n"
            f";磨削余量: {tG}mm;\n"
            f";Y轴加工位置: {Ypos:.4f};\n"
            f";****************\n"
            f";软件版本:1.0.0\n"
            f";生成日期:{current_datetime_str}\n"
            f";****************\n"
            "DEF REAL currentValue;\n"
            "STOPRE;\n"
            "currentValue = R210;\n\n"
        )
        f.write(header)

        # 条件跳转部分
        for i, val in enumerate(step_values):
            label = "CV_" + f"{val:.2f}".replace(".", "_")
            if i == 0:
                # 第一条条件
                line = f"IF (currentValue>={val:.2f}) GOTOF {label};"
            else:
                prev_val = step_values[i - 1]
                line = f"IF (currentValue<{prev_val:.2f}) AND (currentValue>={val:.2f}) GOTOF {label};"
            f.write(line + "\n")
        # 最后一条
        f.write("IF (currentValue<0) GOTOF CV_0_00;\n\n")

        # 分别写入每个标签块
        for val in step_values:
            label = "CV_" + f"{val:.2f}".replace(".", "_")

            # 计算当前 R1
            current_R1 = R1 - val
            current_X0 = X0

            (start_x, start_z, start_b, start_c,
             end_x, end_z,   end_b,   end_c,
             center_x, center_z, center_r) = calculate_grinding_parameters(
                current_R1, alfa1, N1, L1, R2, H1,
                pA1, hN, L21, L22, current_X0, Z0, B0, C0
            )

            block = (
                f"{label}:\n"
                ";基准圆心坐标\n"
                f"R221={center_x:10.4f}        ;X坐标\n"
                f"R222={center_z:10.4f}        ;Z坐标\n"
                f"R223={center_r:10.4f}        ;圆弧半径\n"
                ";============\n"
                ";磨削起点坐标\n"
                f"R231={start_x:10.4f}        ;X坐标\n"
                f"R232={start_z:10.4f}        ;Z坐标\n"
                f"R233={start_b:10.4f}        ;B绝对角度\n"
                f"R234={start_c:10.4f}        ;C绝对角度\n"
                ";============\n"
                ";磨削终点坐标\n"
                f"R241={end_x:10.4f}        ;X坐标\n"
                f"R242={end_z:10.4f}        ;Z坐标\n"
                f"R243={end_b:10.4f}        ;B绝对角度\n"
                f"R244={end_c:10.4f}        ;C轴终端增量角度\n"
                "RET\n\n"
            )
            f.write(block)

    print(f"输出已完成，请查看 {filename} 文件。")

def open_directory(path):
    if platform.system() == "Windows":
        os.startfile(path)
    elif platform.system() == "Darwin":
        subprocess.Popen(["open", path])
    else:
        subprocess.Popen(["xdg-open", path])

def main():
    """
    主函数：
      1. 用户输入砂轮直径范围和步距
      2. 针对每个直径，计算并输出对应的文本文件
      3. 所有文件统一放在文件夹下
    """
    # ====== 通用参数定义（可根据实际情况修改） ======
    R1  = 110.05            # 涡轮分度圆半径
    alfa1 = 21 * 2        # 涡轮加工角度（度）
    N1  = 60                # 涡轮齿数
    L1  = 606.23            # B轴中心距离砂轮中心的 X 向距离
    H1  = 3.40              # 砂轮外沿距离分度圆点的高度
    pA1 = 7.0               # 蜗杆螺旋角度
    hN  = 1                 # 蜗杆旋向 (1右旋 -1左旋)
    L21 = 75.1328           # A轴中心距离砂轮中心的Z向距离
    L22 = -1.145            # B轴中心距离A轴中心
    a0Y = -0.0124          # A轴0度时Y轴和工件中心基准坐标

    # ======== 磨削中点基准坐标和角度 =========
    X0  = 0                 # 磨削中点分度圆点2砂轮对刀点X坐标
    Z0  = 0                 # 磨削中点分度圆点2砂轮对刀点Z坐标
    B0  = 0                 # 磨削中点B轴角度
    C0  = 168.01                 # 磨削中点C轴角度

    # ====== 生成数据参数 ===================
    tG  = 1.2               # 磨削总深度
    fR  = 0.02              # 磨削步进（轴向余量步进）

    min_diameter = 555.0    # 最小砂轮直径
    max_diameter = 557.0    # 最大砂轮直径
    step_diameter = 0.1     # 砂轮步距

    # ====================================

    R1_str = f"{R1:.4f}".replace('.', '_')

    # 确保输出文件夹存在
    output_folder = f"C:/Users/Marco Nie/Downloads/PC_{R1_str}"

    # 删除文件夹
    if os.path.exists(output_folder) and os.path.isdir(output_folder):
        shutil.rmtree(output_folder)

    # 新建文件夹
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Y 轴加工时的坐标
    Ypos = a0Y - math.sin(math.radians(pA1)) * L21 * hN

    # 准备参数字典（先放部分通用）
    param_dict = {
        "R1": R1, "alfa1": alfa1, "N1": N1,
        "L1": L1,             "H1": H1,
        "pA1": pA1,           "hN": hN,
        "L21": L21,           "L22": L22,
        "X0": X0,             "Z0": Z0,
        "B0": B0,             "C0": C0,
        "tG": tG,             "fR": fR,
        "Ypos": Ypos
    }

    # 在指定直径范围内循环
    current_d = min_diameter
    while current_d <= max_diameter + 1e-9:
        # 当前砂轮半径
        R2 = current_d / 2.0
        param_dict["R2"] = R2

        # 生成该直径下的步进序列
        step_values = generate_step_values(tG, fR)

        # 输出文件名，可根据需求自定义
        current_d_str = f"{current_d:.4f}".replace('.', '_')
        output_filename = os.path.join(
            output_folder, f"DIA_{current_d_str}.SPF"
        )

        # 写入输出文件
        write_output_file(output_filename, step_values, param_dict)

        # 下一个直径
        current_d += step_diameter

    print("所有砂轮直径对应的文本文件均已生成。")
    open_directory(output_folder)

if __name__ == "__main__":
    main()
