#!/bin/bash

if [ -z $1 ]; then
  exit 1
fi

FIELD=$1


rm -f frame_????.png
cnt=0

FIELD_NAME=$FIELD

if [ "$FIELD" == 'Temp' ]; then
  FIELD_NAME='Температура'
fi

if [ "$FIELD" == 'Rain' ]; then
  FIELD_NAME='Превалявания'
fi


if [ "$FIELD" == 'Snow' ]; then
  FIELD_NAME='Снеговалеж'
fi


if [ "$FIELD" == 'Cloud' ]; then
  FIELD_NAME='Облачна покривка'
fi

if [ "$FIELD" == 'Humid' ]; then
  FIELD_NAME='Мъгла'
fi

for file in `ls -1 ${FIELD}_????.png`
do
  printf 'Processing %s ... ' $file
  convert $file  -resize x720 -background gray -gravity center -extent 1280x720 -font DejaVu-Sans-Bold -pointsize 30  -encoding Unicode -gravity northwest -fill black -annotate 0 "${FIELD_NAME}" `printf 'frame_%04d.png' $cnt`
  printf 'Done\n'
  let cnt=cnt+1
done

if [ $? -eq 0 ]; then
  ffmpeg -y -r 2/1 -i frame_%04d.png -c:v libx264 -preset slow -pix_fmt yuv420p -refs:v 4 -qp 4 -r 25 $FIELD.mp4
  ffmpeg -i $FIELD.mp4 -q:v 10 -c:v libtheora -c:a libvorbis $FIELD.ogv
  ffmpeg -i $FIELD.mp4 -q:v 10 -c:v libvpx -c:a libvorbis $FIELD.webm
fi

