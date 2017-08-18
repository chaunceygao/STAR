# STAR
STAR code for iLIDS-VID dataset

STAR: Superpixel-based Temporally Aligned Representation for Video-based Person Re-identification

Code to accompany the paper:
	Superpixel-based Temporally Aligned Representation for Video-based Person Re-identification (submitted)

Contact: Changxin Gao <cgao@hust.edu.cn>

# Usage
1 Download the code and the models of Deep Decompositional Network for person segmentation, from http://mmlab.ie.cuhk.edu.hk/projects/luoWTiccv2013DDN/index.html, and unzip it to '.\utils\PedParsing'

2 Download the iLIDS-VID dataset from http://www.eecs.qmul.ac.uk/~xiatian/downloads_qmul_iLIDS-VID_ReID_dataset.html, and modify the path in 'STAR_parameterInitial.m'

3 run 'STAR.m'.

# Note
This code is setting for iLIDS-VID dataset, if you want to evaluate it on other datsets, please change the setting in the file 'STAR_parameterInitial.m'

To easily reproduce our results, we also give our intermediate and final results. They can be downloaded from http://pan.baidu.com/s/1kUNiZDt.

If you have any questions, please contact me: cgao@hust.edu.cn
