msg=`date`
msg="routine backup on $msg"

git add .
git commit -m "$msg"
git push -u origin main
~ 
