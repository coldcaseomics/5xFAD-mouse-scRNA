## 개요
마우스 5xFAD 알츠하이머 모델의 단일핵 RNA 시퀀싱 데이터 분석 파이프라인

- 데이터: GSE140510 (GEO)
- 비교: WT (정상) vs WT_5XFAD (알츠하이머 모델)
- 분석 도구: R / Seurat v5

## 데이터 다운로드
1. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140510 접속
2. GSE140510_RAW.tar 다운로드
3. ./data/ 폴더에 압축 해제

## 폴더 구조
5xFAD_project/
├── 5xFAD_snRNAseq_WT_vs_5xFAD_simple.R
├── README.md
├── plots/
└── rds/


## 분석 단계
1. 데이터 로드 및 QC (nFeature 200~6000, nCount > 500, percent.mt < 5%)
2. SCTransform 정규화
3. PCA → UMAP (PC 30개 사용)
4. 클러스터링 (resolution = 0.5)
5. FindAllMarkers로 세포 타입 annotation
6. WT vs 5xFAD 세포 비율 비교
7. DEG 분석 (미세아교세포, 별아교세포)

## 주요 결과
- 미세아교세포: WT 3.96% → 5xFAD 7.62% (2배 증가)
- 별아교세포: WT 6.89% → 5xFAD 2.05% (3배 감소)
- 미세아교세포 DEG: Cst7, Lpl, Igf1 과발현 (DAM 마커)

## 실행 방법
1. setwd() 경로를 본인 경로로 수정
2. 스크립트 순서대로 실행