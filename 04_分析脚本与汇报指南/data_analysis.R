# ==============================================================================
# 双向孟德尔随机化分析脚本 (Bidirectional Mendelian Randomization)
# 暴露与结局：载脂蛋白B (ApoB) <-> 膝关节骨关节炎 (Knee Osteoarthritis)
# 核心技术：基于 .tbi 索引局部读取与 OpenGWAS API 在线抓取以优化内存分配
# ==============================================================================
library(VariantAnnotation)
library(gwasvcf)
library(gwasglue)
library(TwoSampleMR)
library(MRPRESSO)
library(forestploter)
library(grid)
library(ggplot2)
library(ieugwasr)
library(plinkbinr)

rm(list = ls())
gc() # 清除工作空间变量并释放内存

# ---------------------------------------------------------
# 1. 核心数据文件与参数配置
# 注意：请确保对应 GWAS 数据的 .tbi 索引文件存在于同一工作目录下
# ---------------------------------------------------------
file_ApoB = "met-d-ApoB.vcf.gz"         
file_KneeOA = "ebi-a-GCST007090.vcf.gz" 

sample_size_ApoB = 115000  
sample_size_KneeOA = 403124 

# ---------------------------------------------------------
# 2. 自定义函数：自动化生成绘图与统计结果
# 功能：生成散点图、漏斗图、留一法图、森林图，并计算及格式化 OR 值
# ---------------------------------------------------------
generate_plots <- function(res, dat, prefix, title_desc) {
  p_scatter = mr_scatter_plot(res, dat)[[1]]
  ggsave(paste0(prefix, "1_Scatter.png"), p_scatter, width=8, height=6, dpi=300)
  
  p_funnel = mr_funnel_plot(mr_singlesnp(dat))[[1]]
  ggsave(paste0(prefix, "2_Funnel.png"), p_funnel, width=8, height=6, dpi=300)
  
  p_loo = mr_leaveoneout_plot(mr_leaveoneout(dat))[[1]]
  ggsave(paste0(prefix, "3_LeaveOneOut.png"), p_loo, width=8, height=8, dpi=300)
  
  OR = generate_odds_ratios(res)
  OR_plot = OR[, c("id.exposure", "id.outcome", "nsnp", "method", "or", "pval", "or_lci95", "or_uci95")]
  OR_plot$`OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", OR_plot$or, OR_plot$or_lci95, OR_plot$or_uci95)
  OR_plot$pval = sprintf("%.3f", OR_plot$pval)
  OR_plot$` ` = paste(rep(" ", 60), collapse = " ")
  
  tm = forest_theme(base_size=10, ci_pch=16, ci_col="black", ci_fill="red", refline_lty="dashed")
  p_forest = forest(OR_plot[, c(1:4, 6, 9, 10)], est=OR_plot$or, lower=OR_plot$or_lci95, upper=OR_plot$or_uci95, ci_column=7, ref_line=1, xlim=c(0.8, 1.2), theme=tm)
  p_forest <- edit_plot(p_forest, which="background", row=c(3), gp=gpar(fill="lightsteelblue1"))
  p_forest <- edit_plot(p_forest, row=0, gp=gpar(fontface="bold"), label=title_desc) 
  
  png(paste0(prefix, "4_Summary_Forest.png"), width=13, height=7, units="in", res=300)
  grid.draw(p_forest); dev.off()
  return(OR)
}


# ==============================================================================
# 第一部分：正向孟德尔随机化分析 (暴露: ApoB -> 结局: 膝关节骨关节炎)
# ==============================================================================
print("====== 开始正向分析：ApoB -> 膝关节骨关节炎 ======")

# 1.1 提取暴露变量 (ApoB) 数据并依据 P < 5e-08 进行强效工具变量筛选
exp_vcf_F = VariantAnnotation::readVcf(file_ApoB)
exp_dat_F = gwasvcf_to_TwoSampleMR(vcf = exp_vcf_F, type = "exposure")
exp_dat_F = subset(exp_dat_F, pval.exposure < 5e-08) 

# 1.2 本地连锁不平衡 (LD) 聚类去重 (Clumping)
exp_clump_F = ieugwasr::ld_clump(dplyr::tibble(rsid=exp_dat_F$SNP, pval=exp_dat_F$pval.exposure), clump_kb=10000, clump_r2=0.001, bfile="1kg.v3/EUR", plink_bin=plinkbinr::get_plink_exe())
exp_clump_F = exp_dat_F[exp_dat_F$SNP %in% exp_clump_F$rsid, ]

# 1.3 计算暴露变量的 F 统计量以评估工具变量强度
exp_clump_F$samplesize.exposure = sample_size_ApoB
R2_F = (exp_clump_F$beta.exposure^2) / ((exp_clump_F$beta.exposure^2) + (exp_clump_F$samplesize.exposure * (exp_clump_F$se.exposure^2)))
exp_clump_F$Fz = R2_F * (exp_clump_F$samplesize.exposure - 2) / (1 - R2_F)

# 1.4 提取结局变量 (膝关节骨关节炎) 数据
# 采用 .tbi 索引按需提取相关 SNP，显著降低内存占用
print("正在通过 .tbi 索引秒级提取结局 SNP，内存占用 < 10MB ...")
out_vcf_F = gwasvcf::query_gwas(vcf = file_KneeOA, rsid = exp_clump_F$SNP)
out_dat_F = gwasvcf_to_TwoSampleMR(vcf = out_vcf_F, type = "outcome")

# 1.5 数据协调 (Harmonization) 并统一等位基因方向
exp_clump_F$id.exposure = "ApoB"
out_dat_F$id.outcome = "Knee OA"
dat_F = harmonise_data(exposure_dat = exp_clump_F, outcome_dat = out_dat_F, action = 2)

# 1.6 执行 MR 分析及结果可视化与导出
res_F = mr(dat_F)
OR_F = generate_plots(res_F, dat_F, "Forward_", "正向：ApoB 对 膝关节骨关节炎 的因果影响")
write.csv(OR_F, "Forward_MR_Result.csv", row.names=F)
print("✅ 正向分析完成！")


# ==============================================================================
# 第二部分：反向孟德尔随机化分析 (暴露: 膝关节骨关节炎 -> 结局: ApoB)
# ==============================================================================
# ---------------------------------------------------------
# 2.1 配置 OpenGWAS API 访问令牌 (Token)
# 用于授权在线获取 GWAS 数据
# ---------------------------------------------------------
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyODE5NDQ3ODU0QHFxLmNvbSIsImlhdCI6MTc3NzcyNjg5OCwiZXhwIjoxNzc4OTM2NDk4fQ.qrkuBj5-en6vrkWCzfbH1ADv7uapulp84dgUwsMlw9WKl46074hkFiY_NpUoqFWcXdp9W4QTMZEO7WbR8MPwXRvff42T7qePSormtg2Gp-XhjWGS6kqee5gWRfnsV5YvUgmunk2-o0bNXMbgkOyD3rCN8rg8PpCcphMWPKCH2BDIVVyhcqPe9nrH8ww6_I__zXTw1ZA9WyBdrEg-d_Mq0cOOiDqgWmvJAoGbFfV6Oo-En0uImzxDTSSjR9wSn2jzbhijP3GkWQcGwnbo5nJINMLGr7VGLLlUcqyFkF1TvIx6pY-jD3QKsQ-HcAFb_FpwzkxIlYEzYEL92XoLerhoXA") 

# ---------------------------------------------------------
# 2.2 执行反向 MR 分析流程
# ---------------------------------------------------------
print("====== 开始反向分析：膝关节骨关节炎 -> ApoB ======")
print("通行证已加载，正在通过在线 API 极速抓取强效 SNP...")

# 2.2.1 提取反向暴露变量 (膝关节骨关节炎)
# 通过在线 API 直接抓取全基因组显著的 SNP，避免加载超大本地文件
print("正在通过在线 API 直接抓取骨关节炎的强效 SNP，保护内存...")
exp_dat_R = extract_instruments(outcomes = "ebi-a-GCST007090", p1 = 5e-08, clump = FALSE)

# 2.2.2 离线 LD 聚类去重 (Clumping)
exp_clump_R = ieugwasr::ld_clump(dplyr::tibble(rsid=exp_dat_R$SNP, pval=exp_dat_R$pval.exposure), clump_kb=10000, clump_r2=0.001, bfile="1kg.v3/EUR", plink_bin=plinkbinr::get_plink_exe())
exp_clump_R = exp_dat_R[exp_dat_R$SNP %in% exp_clump_R$rsid, ]

# 2.2.3 计算反向暴露变量的 F 统计量
exp_clump_R$samplesize.exposure = sample_size_KneeOA
R2_R = (exp_clump_R$beta.exposure^2) / ((exp_clump_R$beta.exposure^2) + (exp_clump_R$samplesize.exposure * (exp_clump_R$se.exposure^2)))
exp_clump_R$Fz = R2_R * (exp_clump_R$samplesize.exposure - 2) / (1 - R2_R)

# 2.2.4 提取反向结局变量 (ApoB) 数据
# 复用 .tbi 索引局部读取技术
print("正在通过 .tbi 索引提取 ApoB 数据...")
out_vcf_R = gwasvcf::query_gwas(vcf = file_ApoB, rsid = exp_clump_R$SNP)
out_dat_R = gwasvcf_to_TwoSampleMR(vcf = out_vcf_R, type = "outcome")

# 2.2.5 反向数据协调 (Harmonization)
exp_clump_R$id.exposure = "Knee OA"
out_dat_R$id.outcome = "ApoB"
dat_R = harmonise_data(exposure_dat = exp_clump_R, outcome_dat = out_dat_R, action = 2)

# 2.2.6 执行反向 MR 分析及结果生成
res_R = mr(dat_R)
OR_R = generate_plots(res_R, dat_R, "Reverse_", "反向：膝关节骨关节炎 对 ApoB 的因果影响")
write.csv(OR_R, "Reverse_MR_Result.csv", row.names=F)
print("✅ 反向分析完成！")

View(OR_F, title="正向结果 (ApoB致病)")
View(OR_R, title="反向结果 (骨关节炎致病)")
print("🎉 双向分析全部圆满结束！")


# ==============================================================================
# 第三部分：质量控制 (QC) - 异质性与水平多效性检验
# ==============================================================================
print("====== 开始生成并保存质检表格 ======")

# 3.1 导出正向分析的异质性与水平多效性检验结果
write.csv(mr_heterogeneity(dat_F), "Forward_Heterogeneity.csv", row.names=F)
write.csv(mr_pleiotropy_test(dat_F), "Forward_Pleiotropy.csv", row.names=F)

# 3.2 导出反向分析的异质性与水平多效性检验结果
write.csv(mr_heterogeneity(dat_R), "Reverse_Heterogeneity.csv", row.names=F)
write.csv(mr_pleiotropy_test(dat_R), "Reverse_Pleiotropy.csv", row.names=F)

print("✅ 所有质检 CSV 表格已成功保存到当前文件夹！")


# ==============================================================================
# 第四部分：质检 (QC) 汇总数据可视化 (生成 Dashboard)
# ==============================================================================
# 若未安装 ggplot2 与 patchwork，请取消下方注释进行安装：
# install.packages("ggplot2")
# install.packages("patchwork")

library(ggplot2)
library(dplyr)
library(patchwork) # 用于多图表拼接排版

print("正在读取本地质检表格...")

# 4.1 读取并处理多效性 (Pleiotropy) 数据，计算 95% 置信区间
f_plei <- read.csv("Forward_Pleiotropy.csv")
r_plei <- read.csv("Reverse_Pleiotropy.csv")

plei_dat <- data.frame(
  Analysis = factor(c("Reverse (Knee OA -> ApoB)", "Forward (ApoB -> Knee OA)"), 
                    levels = c("Reverse (Knee OA -> ApoB)", "Forward (ApoB -> Knee OA)")),
  Intercept = c(r_plei$egger_intercept[1], f_plei$egger_intercept[1]),
  SE = c(r_plei$se[1], f_plei$se[1]),
  Pval = c(r_plei$pval[1], f_plei$pval[1])
)
plei_dat$LCI <- plei_dat$Intercept - 1.96 * plei_dat$SE
plei_dat$UCI <- plei_dat$Intercept + 1.96 * plei_dat$SE

# 4.2 读取并处理异质性 (Heterogeneity) 数据，提取 IVW 方法 P 值
f_het <- read.csv("Forward_Heterogeneity.csv")
r_het <- read.csv("Reverse_Heterogeneity.csv")

f_het_ivw <- f_het %>% filter(method == "Inverse variance weighted") %>% pull(Q_pval)
r_het_ivw <- r_het %>% filter(method == "Inverse variance weighted") %>% pull(Q_pval)

het_dat <- data.frame(
  Analysis = factor(c("Reverse (Knee OA -> ApoB)", "Forward (ApoB -> Knee OA)"),
                    levels = c("Reverse (Knee OA -> ApoB)", "Forward (ApoB -> Knee OA)")),
  Pval = c(r_het_ivw, f_het_ivw)
)
het_dat$Significance <- ifelse(het_dat$Pval < 0.05, "Heterogeneity (P<0.05)", "No Heterogeneity (P>0.05)")

# 4.3 绘制左侧子图：水平多效性森林图
p_plei <- ggplot(plei_dat, aes(x = Intercept, y = Analysis)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d62728", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2, color = "#2ca02c", linewidth = 1.2) +
  geom_point(color = "#2ca02c", size = 4) +
  geom_text(aes(label = paste0("P = ", round(Pval, 3))), vjust = -1.5, color = "#2ca02c", fontface = "bold") +
  theme_minimal(base_size = 14) +
  labs(title = "Horizontal Pleiotropy Test",
       subtitle = "MR-Egger Intercept (Crossing 0 = No Pleiotropy)",
       x = "Intercept (95% CI)", y = "") +
  theme(plot.title = element_text(face = "bold"))

# 4.4 绘制右侧子图：异质性评估柱状图
p_het <- ggplot(het_dat, aes(x = Analysis, y = Pval, fill = Significance)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.85) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#d62728", linewidth = 1) +
  geom_text(aes(label = sprintf("P = %.2e", Pval)), vjust = -0.8, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Heterogeneity (P<0.05)" = "#ff7f0e", "No Heterogeneity (P>0.05)" = "#1f77b4")) +
  theme_minimal(base_size = 14) +
  labs(title = "Heterogeneity Test",
       subtitle = "IVW method (Dashed line at P=0.05)",
       x = "", y = "P-value") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  coord_flip() # 翻转坐标轴以对齐整体视觉结构

# 4.5 拼接图表并导出为高分辨率 PNG 及矢量化 PDF 文件
final_plot <- p_plei + p_het + 
  plot_annotation(title = 'Mendelian Randomization Quality Control (QC) Summary',
                  theme = theme(plot.title = element_text(size = 18, face = 'bold', hjust = 0.5)))

ggsave("Final_QC_Dashboard.png", final_plot, width = 14, height = 6, dpi = 300, bg = "white")
ggsave("Final_QC_Dashboard.pdf", final_plot, width = 14, height = 6)

print("✅ 绘图大功告成！请在文件夹查看 Final_QC_Dashboard.png 和 .pdf 文件！")


# ==============================================================================
# 第五部分：输出最终的协调后数据 (Harmonized Data / SNP 详情明细)
# 包含对齐后的 Allele、Beta、SE、P 值及 F 统计量，用于论文的 Supplementary Tables
# ==============================================================================

# 5.1 导出正向分析的 SNP 特征参数表 (ApoB -> 膝关节骨关节炎)
write.csv(dat_F, "Forward_SNP_Details.csv", row.names = FALSE)

# 5.2 导出反向分析的 SNP 特征参数表 (膝关节骨关节炎 -> ApoB)
write.csv(dat_R, "Reverse_SNP_Details.csv", row.names = FALSE)