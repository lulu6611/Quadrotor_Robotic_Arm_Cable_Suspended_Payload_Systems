
# 四旋翼–机械臂–缆绳负载系统的物理建模与数据驱动混合控制（论文配套仓库）

# Physical Modeling and Data-Driven Hybrid Control for Quadrotor–Robotic-Arm Cable-Suspended Payload Systems (Project Repository)

---

## 简介

本仓库用于展示论文《Physical Modeling and Data-Driven Hybrid Control for Quadrotor–Robotic-Arm Cable-Suspended Payload Systems》中涉及的部分示例代码、MuJoCo 模型、轨迹生成脚本、控制器实现、实验数据样例和仿真视频。
The repository provides example code, MuJoCo models, trajectory generators, controller implementations, experimental data samples, and simulation videos associated with the paper *“Physical Modeling and Data-Driven Hybrid Control for Quadrotor–Robotic-Arm Cable-Suspended Payload Systems.”*

仓库内容旨在帮助研究者复现论文主要实验结果，并为进一步研究提供可扩展的基础框架。
The materials are intended to help researchers reproduce the main experimental results and provide an extensible framework for further development.

本仓库中的代码仅用于学术参考，不构成可直接部署的飞控系统。
All code is for academic use only and is not intended as deployable flight-control firmware.

<div align="center">
  <img src="system structure drawing.png" alt="UAV结构图" width="700">
  <p><em>图: 四旋翼-机械臂-悬挂负载系统结构图</em></p>
  <p><em>Figure: Schematic diagram of the quadrotor-robotic-arm-suspended load system</em></p>
</div>

---

## 仓库结构说明

##  📁 Repository Structure

### 📁 docs

用于存放论文 PDF、补充说明文档以及参考文献 BibTeX 文件。
Contains the paper PDF, supplementary notes, and bibliographic references.

---

### 📁 model

包含四旋翼、两段机械臂、缆绳与负载系统的 MuJoCo 模型，包括：
Includes MuJoCo models of the quadrotor, two-link robotic arm, cable, and payload:

* uav_MW.xml 主模型文件
  uav_MW.xml as the main MuJoCo simulation model
* 📁 meshes 文件夹：机体、机械臂、传感器、负载等三维网格
  meshes folder containing 3D geometry of the UAV body, arm links, sensors, and payload

---

### 📁 src

核心代码目录，分为以下子模块：
The main source code directory, organized into the following components:

#### 不同绳长、质量、轨迹下的飞行原始数据  📁Identification-ready-raw-data

Raw flight and simulation data across different cable lengths, payload masses, and trajectories.

#### 数据驱动的辨识脚本

Data-driven identification scripts for model parameter estimation.

#### 不同实验的测试轨迹

Test trajectories used across experiments.

---

### 📦 quadrotor-arm-cable-hybrid-control.zip

直接运行的仿真程序，轨迹可以自己设置，已包含辨识结果。
A ready-to-run simulation environment that allows custom trajectory inputs and includes identified model parameters.

---

### 📁 videos

存放对应实验的仿真视频，包括：
Contains simulation videos for each experiment, including:

* 阶跃响应视频（MuJoCo 录制 + MATLAB 录制）
  Step-response videos (MuJoCo capture + MATLAB visualization)
* 八字轨迹视频（MuJoCo 录制 + MATLAB 录制）
  Figure-eight trajectory videos (MuJoCo + MATLAB)
* 机械臂协同消摆视频（MuJoCo 录制 + MATLAB 录制）
  Cooperative arm swing suppression videos
* 参数鲁棒性测试视频（MuJoCo 录制 + MATLAB 录制）
  Parametric robustness testing videos
* 在线质量估计视频（MATLAB 录制）
  Online payload mass estimation videos

---

## 使用说明

## Instructions for Use

需要提前在电脑上安装 MATLAB 2024b 及以上版本，并且安装支持 MuJoCo 的 MATLAB 工具箱。
MATLAB 2024b or later is required, along with the MATLAB–MuJoCo interface toolbox.

工具箱安装链接：
Toolbox installation link:
[https://ww2.mathworks.cn/matlabcentral/fileexchange/128028-simulink-blockset-for-mujoco-simulator/](https://ww2.mathworks.cn/matlabcentral/fileexchange/128028-simulink-blockset-for-mujoco-simulator/)

---

### 克隆仓库

### Clone the Repository

```bash
git clone https://github.com/lulu6611/Quadrotor_Robotic_Arm_Cable_Suspended_Payload_Systems.git
cd quadrotor-arm-cable-hybrid-control
```

---

### 修改模型读取路径

### Configure Model File Paths

0. 将仓库中的📦 quadrotor-arm-cable-hybrid-control.zip 保存在电脑上
   Save 📦 quadrotor-arm-cable-hybrid-control.zip locally on your computer.

1. 文件夹内容应如下，确保 uav_MW.xml 与 📁meshes 文件夹放在同一目录下：
   Ensure the folder structure matches the list below and uav_MW.xml is placed together with the meshes folder.

```
quadrotor-arm-cable-hybrid-control/
├── uav_MW.xml
├── Uav_YesArm_MCG_1st_control_192.slx
├── meshes/
│   ├── army.STL
│   ├── armz.STL
│   ├── base_link.STL
│   ├── p1.STL
│   ├── p2.STL
│   ├── p3.STL
│   ├── p4.STL
│   ├── winch.STL
│   ├── d435i_0.obj
│   ├── d435i_1.obj
│   ├── d435i_2.obj
│   ├── d435i_3.obj
│   ├── d435i_4.obj
│   ├── d435i_5.obj
│   ├── d435i_6.obj
│   ├── d435i_7.obj
│   └── d435i_8.obj
```

2. 打开 Uav_YesArm_MCG_1st_control_18.slx，将其中“MuJoCo interaction module”模块的 XML 文件路径设置为 uav_MW.xml
   Open Uav_YesArm_MCG_1st_control_18.slx and configure the XML file path inside the “MuJoCo interaction module” block to point to uav_MW.xml.

3. 运行示例仿真 Uav_YesArm_MCG_1st_control_18.slx
   Run the demonstration model.

---

## 引用

## Citation

如果您在研究中使用本仓库，请引用论文：
If you use this repository in your research, please cite:

```bibtex
@article{Lu2025QuadrotorArmCableHybrid,
  author  = {Lu, Lu and Xiao, Qihua and Zhou, Shikang and Wang, Xinhai and Meng, Yunhe},
  title   = {Physical Modeling and Data-Driven Hybrid Control for Quadrotor--Robotic-Arm Cable-Suspended Payload Systems},
  journal = {Drones},
  year    = {2025},
  volume  = {X},
  number  = {X},
  pages   = {X--X},
  doi     = {10.3390/drones10XXXX}
}
```

---

## 联系方式

## Contact

如有问题或希望合作，可通过邮件联系：
For questions or collaboration inquiries, please contact:

Lu Lu: [lulu28@mail2.sysu.edu.cn](mailto:lulu28@mail2.sysu.edu.cn) or
       [lulu28@mail3.sysu.edu.cn](mailto:lulu28@mail3.sysu.edu.cn)

