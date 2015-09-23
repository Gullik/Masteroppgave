TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    ../src/*.c \
    ../src/*.h \
    ../src/*.cpp \
    ../src/accel.c \
    ../src/collisions.c \
    ../src/collisions_constant.c \
    ../src/diagn.c \
    ../src/dustg.c \
    ../src/flux.c \
    ../src/gauss.c \
    ../src/generate.c \
    ../src/grid.c \
    ../src/input.c \
    ../src/main.c \
    ../src/photons.c \
    ../src/restart.c \
    ../src/shortcuts.c \
    ../src/spherical.c \
    ../src/fmg/fmg.c \
    ../src/fmg/fmg_P.c \
    ../src/fmg/nrutil.c

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    ../src/const.h \
    ../src/funct.h \
    ../src/fmg/nrutil.h

