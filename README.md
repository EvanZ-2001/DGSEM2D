# DGSEM2D

V2版本的代码中squeeze函数占用了太多的运行时间，V3版本（当前版本）主要对调用squeeze的部分进行改变，将单元序号改为最后一项索引，避免使用squeeze，对于仍然需要squeeze的地方，改用reshape函数

