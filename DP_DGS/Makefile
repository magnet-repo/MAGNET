CC = gcc
DEBUG = FALSE

CFLAGS = -O3 -march=native -fomit-frame-pointer 
WARNINGS = -Wall -Wextra -Wimplicit-function-declaration -Wredundant-decls -Wmissing-prototypes -Wstrict-prototypes -Wundef -Wno-deprecated-declarations
LDFLAGS = -lm

ifeq ($(DEBUG),TRUE)
  CFLAGS += -DDEBUG
endif

SRC   := build.c
OBJDIR  := objs
OBJ     := $(SRC:%.c=$(OBJDIR)/%.o)
TARGET := DP-DGS

all: clean $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(WARNINGS) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: %.c | $(OBJDIR)
	$(CC) $(CFLAGS) $(WARNINGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(TARGET)
