
cd ~/project1/github_clean


git init


git lfs install


git lfs track "*.zip"


git add .gitattributes
git commit -m "Track .zip files with Git LFS"


git add --all


git status


git commit -m "Add all project files"


git remote add origin https://github.com/LinResearchGroup-NCSU/IDEA_Model.git


git push -u origin main

