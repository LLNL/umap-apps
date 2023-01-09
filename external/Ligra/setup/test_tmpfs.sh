#!/bin/bash


#function setuptmpfs {
  if [ ! -d /mnt/tmpfs ]; then
    mkdir -p /mnt/tmpfs
  fi

  # Unmount / Reset of already mounted
  fs=`stat -f -c '%T' /mnt/tmpfs`

  if [ "$fs" = "tmpfs" ]; then
    echo "Resetting tmpfs"
    umount /mnt/tmpfs
  fi

  fs=`stat -f -c '%T' /mnt/tmpfs`
  if [ "$fs" != "tmpfs" ]; then
    if [ ! -d /mnt/tmpfs ]; then
      mkdir -p /mnt/tmpfs
    fi
    chmod go+rwx /mnt/tmpfs
    mount -t tmpfs -o size=48g tmpfs /mnt/tmpfs
    fs=`stat -f -c '%T' /mnt/tmpfs`
    echo "/mnt/tmpfs mounted as: $fs"
  else
    echo "Unable to reset /mnt/tmpfs, exiting"
    exit 1
  fi
#}
