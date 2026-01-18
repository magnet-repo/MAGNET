# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.20.0)

find_package(Zephyr REQUIRED HINTS $ENV{ZEPHYR_BASE})
project(MAGNET_m4)

enable_language(ASM)

file(GLOB sources "src/*.c" "src/*.s" "src/*.S")

target_sources(app PRIVATE src/main.c ${sources})
