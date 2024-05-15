# 导航到你的项目目录
cd ~/project1/github_clean

# 初始化 Git 仓库
git init

# 初始化 Git LFS
git lfs install

# 设置 Git LFS 跟踪大文件类型
git lfs track "*.zip"

# 提交 .gitattributes 文件
git add .gitattributes
git commit -m "Track .zip files with Git LFS"

# 添加所有未跟踪的文件和文件夹
git add --all

# 确认所有文件都已添加
git status

# 提交更改
git commit -m "Add all project files"

# 添加远程仓库
git remote add origin https://github.com/LinResearchGroup-NCSU/IDEA_Model.git

# 推送到远程仓库的 main 分支
git push -u origin main

